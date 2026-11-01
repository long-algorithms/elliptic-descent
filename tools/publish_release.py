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


def gh(*args, optional=False, input_data=None):
    p = subprocess.run(['gh', *args], input=input_data, capture_output=True, check=False)
    if p.returncode:
        if optional and b'HTTP 404' in p.stderr:
            return None
        raise RuntimeError(p.stderr.decode(errors='replace'))
    return p.stdout


def api(path, optional=False, *, method=None, payload=None):
    args = ['api', path]
    if method is not None:
        args.extend(('--method', method))
    data = None
    if payload is not None:
        args.extend(('--input', '-'))
        data = json.dumps(payload).encode('utf-8')
    result = gh(*args, optional=optional, input_data=data)
    return None if result is None else json.loads(result)


def find_release(repo, tag):
    """Find published OR draft releases, retaining their numeric identity.

    GitHub's release-by-tag endpoint does not reliably return unpublished
    releases. Authenticated listing includes them. Exhaust the list to reject
    ambiguous draft tags, rather than accidentally creating another draft.
    """
    found, page = [], 1
    while True:
        batch = api(f'repos/{repo}/releases?per_page=100&page={page}')
        require(isinstance(batch, list), 'invalid release listing')
        found.extend(r for r in batch if r.get('tag_name') == tag)
        if len(batch) < 100:
            break
        page += 1
    require(len(found) <= 1, 'multiple releases have the requested tag')
    return found[0] if found else None


def release_identity(release, tag, sha, *, draft):
    require(type(release.get('id')) is int and release['id'] > 0
            and release.get('tag_name') == tag
            and release.get('target_commitish') == sha
            and release.get('draft') is draft
            and release.get('prerelease') is False,
            'release identity or publication state changed')
    if draft:
        require(release.get('published_at') is None, 'draft was previously published')


def replacement_allowed(repo, release, tag, sha, replacement):
    """Authorize only an explicitly named, unpublished, untagged ancestor draft.

    This does not mutate anything and does not waive any new-revision audit.
    Both the draft ID and former exact revision must be explicitly supplied.
    """
    match = re.fullmatch(r'([1-9][0-9]*)@([0-9a-f]{40})', replacement or '')
    require(match is not None, 'existing draft targets another revision; explicit recovery required')
    identifier, old_sha = int(match[1]), match[2]
    release_identity(release, tag, old_sha, draft=True)
    require(release['id'] == identifier and old_sha != sha, 'wrong recovery draft')
    require(tag_commit(repo, tag) is None, 'refusing to retarget an existing release tag')
    comparison = api(f'repos/{repo}/compare/{old_sha}...{sha}')
    require(comparison.get('status') == 'ahead'
            and comparison.get('merge_base_commit', {}).get('sha') == old_sha,
            'replacement revision is not a descendant of the audited draft')


def check_uploaded_assets(release, output):
    files = {p.name: p for p in output.iterdir() if p.is_file()}
    assets = release.get('assets', [])
    require(len(assets) == len(files) and {a['name'] for a in assets} == set(files),
            'incomplete or duplicate asset upload')
    for asset in assets:
        content = files[asset['name']].read_bytes()
        require(asset.get('state') == 'uploaded' and asset.get('size') == len(content)
                and asset.get('digest') == 'sha256:'+hashlib.sha256(content).hexdigest(),
                'uploaded asset checksum mismatch')


def stage_and_publish(repo, sha, tag, output, existing, replacement=None):
    """Called only after all exact-head workflow, evidence and package audits.

    Work by numeric release ID. A recovery can only retarget an explicitly
    named unpublished ancestor draft, after rechecking it and main. No public
    release, existing tag or unknown draft is overwritten.
    """
    require(api(f'repos/{repo}/git/ref/heads/main')['object']['sha'] == sha,
            'main moved during audit')
    if existing is None:
        require(find_release(repo, tag) is None, 'a release appeared during audit')
        gh('release', 'create', tag, '--repo', repo, '--target', sha, '--draft',
           '--title', 'v0.5.0 — explicit 4-Selmer groups over Q',
           '--notes-file', str(ROOT/'docs/release-notes-v0.5.0.md'))
        existing = find_release(repo, tag)
        require(existing is not None, 'created draft was not found')
    else:
        current = api(f'repos/{repo}/releases/{existing["id"]}')
        require(current['id'] == existing['id'], 'draft ID changed')
        if current.get('target_commitish') != sha:
            replacement_allowed(repo, current, tag, sha, replacement)
            current = api(f'repos/{repo}/releases/{current["id"]}', method='PATCH',
                          payload={'target_commitish': sha})
        existing = current
    release_identity(existing, tag, sha, draft=True)
    identifier = existing['id']
    require(tag_commit(repo, tag) in (None, sha), 'release tag changed during audit')
    gh('release', 'upload', tag, '--repo', repo, '--clobber', *map(str, sorted(output.iterdir())))
    draft = api(f'repos/{repo}/releases/{identifier}')
    release_identity(draft, tag, sha, draft=True)
    check_uploaded_assets(draft, output)
    require(api(f'repos/{repo}/git/ref/heads/main')['object']['sha'] == sha,
            'main moved before publication')
    require(tag_commit(repo, tag) in (None, sha), 'release tag changed before publication')
    api(f'repos/{repo}/releases/{identifier}', method='PATCH',
        payload={'draft': False, 'make_latest': 'true'})
    final = api(f'repos/{repo}/releases/{identifier}')
    release_identity(final, tag, sha, draft=False)
    check_uploaded_assets(final, output)
    require(tag_commit(repo, tag) == sha, 'published release tag mismatch')
    return final


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


def publish(repo, sha, *, replacement=None):
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
    existing = find_release(repo, tag)
    if existing is not None:
        if existing.get('target_commitish') != sha:
            replacement_allowed(repo, existing, tag, sha, replacement)
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
        final = stage_and_publish(repo, sha, tag, output, existing, replacement)
        print(json.dumps({'published': True, 'tag': tag, 'revision': sha,
                          'url': final['html_url'], 'audit': report}, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--replace-draft', help='Explicit unpublished draft recovery: numeric ID@old commit SHA')
    args = parser.parse_args()
    publish(args.repository, args.revision, replacement=args.replace_draft)
