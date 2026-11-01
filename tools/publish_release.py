"""Publish only a fully audited, exact main-branch revision to GitHub Releases.

Runs inside the explicit release workflow, using its short-lived GH_TOKEN.
No PyPI upload, secret discovery, tag rewriting, or branch modification occurs.
An incomplete gate makes this invocation a no-op; a failed gate is an error.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipfile

from release_gate import VERSION, audit_directory, require, selected_runs
from check_distribution import check_distributions

ROOT = Path(__file__).resolve().parents[1]


def gh(*args, optional=False):
    p = subprocess.run(['gh', *args], capture_output=True, check=False)
    if p.returncode:
        if optional and b'HTTP 404' in p.stderr:
            return None
        raise RuntimeError(p.stderr.decode(errors='replace'))
    return p.stdout


def api(path, optional=False):
    data = gh('api', path, optional=optional)
    return None if data is None else json.loads(data)


def pages(path, key):
    result, page = [], 1
    while True:
        separator = '&' if '?' in path else '?'
        data = api(f'{path}{separator}per_page=100&page={page}')[key]
        result.extend(data)
        if len(data) < 100:
            return result
        page += 1


def extract_archive(data, destination):
    """Reject traversal, symlinks and duplicate paths in downloaded artifacts."""
    import io
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        seen = set()
        for member in archive.infolist():
            p = Path(member.filename)
            require(not p.is_absolute() and '..' not in p.parts
                    and '\\' not in member.filename and p.as_posix() not in seen,
                    'unsafe/duplicate artifact member')
            require((member.external_attr >> 16) & 0o170000 != 0o120000, 'artifact symlink rejected')
            seen.add(p.as_posix())
        archive.extractall(destination)


def tag_commit(repo, tag):
    ref = api(f'repos/{repo}/git/ref/tags/{tag}', optional=True)
    if ref is None:
        return None
    obj = ref['object']
    for _ in range(8):
        if obj['type'] == 'commit':
            return obj['sha']
        require(obj['type'] == 'tag', 'release tag has an invalid object type')
        obj = api(f'repos/{repo}/git/tags/{obj["sha"]}')['object']
    raise ValueError('excessively nested annotated tag')


def publish(repo, sha):
    require(re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', repo) is not None, 'invalid repository')
    require(re.fullmatch(r'[0-9a-f]{40}', sha) is not None, 'invalid revision')
    require(os.environ.get('GITHUB_REPOSITORY') == repo, 'wrong workflow repository')
    actual = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    require(actual == sha, 'checkout is not the requested revision')
    version = re.search(r'^version = "([^"]+)"', (ROOT/'pyproject.toml').read_text(), re.M).group(1)
    if version != VERSION:
        print('No v0.5.0 release requested by this source revision.')
        return
    require(api(f'repos/{repo}/git/ref/heads/main')['object']['sha'] == sha, 'main moved; re-audit the new revision')
    runs = selected_runs(pages(f'repos/{repo}/actions/runs?head_sha={sha}&event=push', 'workflow_runs'), sha)
    if runs is None:
        print('Release not ready: required exact-head workflows are still missing or running.')
        return
    tag = 'v'+version
    existing_tag = tag_commit(repo, tag)
    require(existing_tag in (None, sha), 'refusing to replace a release tag')
    existing = api(f'repos/{repo}/releases/tags/{tag}', optional=True)
    if existing is not None:
        require(existing.get('target_commitish') == sha, 'existing release targets another revision')
        if not existing['draft']:
            require(existing_tag == sha, 'published release tag mismatch')
            print(f'{tag} is already published at the verified revision {sha}.')
            return

    with tempfile.TemporaryDirectory(prefix='elliptic-release-') as tmp:
        work = Path(tmp)
        artifacts, evidence, output = work/'artifacts', work/'evidence', work/'output'
        for p in (artifacts, evidence, output):
            p.mkdir()
        provenance = {'revision': sha, 'version': version, 'runs': runs, 'artifacts': []}
        # Retain all small exact evidence. The optional third-party accelerator
        # wheel is not a project release asset and need not be redistributed.
        for run in runs.values():
            for a in pages(f'repos/{repo}/actions/runs/{run["id"]}/artifacts', 'artifacts'):
                if a['name'] == 'exact-accelerator-wheel':
                    continue
                require(not a['expired'] and re.fullmatch(r'[A-Za-z0-9_.-]+', a['name']) is not None,
                        'expired/unsafe artifact')
                data = gh('api', f'repos/{repo}/actions/artifacts/{a["id"]}/zip')
                digest = hashlib.sha256(data).hexdigest()
                require(a.get('digest') == 'sha256:'+digest, 'artifact digest mismatch')
                destination = artifacts/a['name']
                require(not destination.exists(), 'duplicate artifact name across workflows')
                extract_archive(data, destination)
                (evidence/(a['name']+'.zip')).write_bytes(data)
                provenance['artifacts'].append({'id': a['id'], 'name': a['name'],
                                                'run': run['id'], 'sha256': digest})
        report = audit_directory(artifacts, ROOT, sha, version)
        manifests = list(artifacts.rglob('release-manifest.json'))
        require(len(manifests) == 1, 'missing/duplicate distribution manifest')
        original = json.loads(manifests[0].read_text())
        require(original['verified'] is True and original['version'] == version
                and original['revision'] == sha, 'distribution provenance mismatch')
        for a in original['artifacts']:
            p = manifests[0].parent/a['file']
            require(p.name == a['file'] and p.is_file() and p.stat().st_size == a['bytes']
                    and hashlib.sha256(p.read_bytes()).hexdigest() == a['sha256'], 'distribution checksum mismatch')
            shutil.copyfile(p, output/p.name)
        os.environ['GITHUB_SHA'] = sha
        checked = check_distributions(ROOT, output)
        require(checked == original, 'fresh distribution audit differs from CI')
        report['distribution'] = checked
        (output/'release-audit.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
        (evidence/'release-audit.json').write_bytes((output/'release-audit.json').read_bytes())
        (evidence/'provenance.json').write_text(json.dumps(provenance, indent=2, sort_keys=True)+'\n')
        subprocess.run(['git', 'archive', '--format=zip', '--prefix=source/',
                        f'--output={evidence / "source.zip"}', sha], cwd=ROOT, check=True)
        bundle = output/f'elliptic-descent-{tag}-verification.zip'
        with zipfile.ZipFile(bundle, 'w', zipfile.ZIP_DEFLATED) as archive:
            for p in sorted(evidence.iterdir()):
                archive.write(p, p.name)
        (output/'SHA256SUMS').write_text(''.join(
            f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n'
            for p in sorted(output.iterdir()) if p.name != 'SHA256SUMS'))
        require(api(f'repos/{repo}/git/ref/heads/main')['object']['sha'] == sha, 'main moved during audit')
        if existing is None:
            gh('release', 'create', tag, '--repo', repo, '--target', sha, '--draft',
               '--title', 'v0.5.0 — explicit 4-Selmer groups over Q',
               '--notes-file', str(ROOT/'docs/release-notes-v0.5.0.md'))
        # Replacing assets is permitted only on this same-commit draft. The
        # release becomes public only after the entire set has uploaded.
        gh('release', 'upload', tag, '--repo', repo, '--clobber', *map(str, sorted(output.iterdir())))
        draft = api(f'repos/{repo}/releases/tags/{tag}')
        require(draft['draft'] and draft['target_commitish'] == sha, 'release changed during upload')
        require({a['name'] for a in draft['assets']} == {p.name for p in output.iterdir()}, 'incomplete asset upload')
        require(api(f'repos/{repo}/git/ref/heads/main')['object']['sha'] == sha, 'main moved before publication')
        gh('release', 'edit', tag, '--repo', repo, '--draft=false', '--latest')
        final = api(f'repos/{repo}/releases/tags/{tag}')
        require(not final['draft'] and not final['prerelease'] and tag_commit(repo, tag) == sha,
                'published release verification failed')
        print(json.dumps({'published': True, 'tag': tag, 'revision': sha,
                          'url': final['html_url'], 'audit': report}, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', required=True)
    parser.add_argument('--revision', required=True)
    args = parser.parse_args()
    publish(args.repository, args.revision)
