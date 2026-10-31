"""Audit built wheel/sdist contents against this source tree and write checksums.

No publishing, credentials or network calls. Run after python -m build.
"""
import argparse
from email.parser import BytesParser
import hashlib
import json
import os
from pathlib import Path
import re
import tarfile
from zipfile import ZipFile


def check_distributions(root, directory):
    project = (root/'pyproject.toml').read_text()
    version = re.search(r'^version = "([^"]+)"', project, re.M).group(1)
    runtime = (root/'src/elliptic_descent/__init__.py').read_text()
    if f'__version__ = "{version}"' not in runtime:
        raise ValueError('runtime and package versions disagree')
    wheels, sdists = list(directory.glob('*.whl')), list(directory.glob('*.tar.gz'))
    if len(wheels) != 1 or len(sdists) != 1:
        raise ValueError('expected exactly one wheel and one source distribution')
    with ZipFile(wheels[0]) as wheel:
        files = wheel.namelist()
        metadata_path = next(name for name in files if name.endswith('.dist-info/METADATA'))
        metadata = BytesParser().parsebytes(wheel.read(metadata_path))
        if metadata['Version'] != version or metadata['Name'] != 'elliptic-descent':
            raise ValueError('wheel metadata differs from project')
        if metadata['Requires-Python'] != '>=3.10' or metadata['License-Expression'] != 'MIT':
            raise ValueError('wheel Python/license metadata is incorrect')
        source_modules = {p.relative_to(root/'src').as_posix(): p for p in (root/'src').rglob('*.py')}
        if {name for name in files if name.endswith('.py')} != set(source_modules):
            raise ValueError('wheel Python modules do not match source modules')
        for name, path in source_modules.items():
            if wheel.read(name) != path.read_bytes():
                raise ValueError(f'wheel module differs: {name}')
    with tarfile.open(sdists[0], 'r:gz') as archive:
        members = {m.name: m for m in archive.getmembers() if m.isfile()}
        prefix = next(iter(members)).split('/')[0]
        required = [root/p for p in ('LICENSE','README.md','CHANGELOG.md','pyproject.toml','MANIFEST.in')]
        for folder, extension in [('src','*.py'),('tests','*.py'),('benchmarks','*.py'),('tools','*.py'),('docs','*.md')]:
            required.extend((root/folder).rglob(extension))
        for path in required:
            name = prefix+'/'+path.relative_to(root).as_posix()
            if name not in members or archive.extractfile(members[name]).read() != path.read_bytes():
                raise ValueError(f'sdist file absent or different: {name}')
    records = []
    for path in (*wheels, *sdists):
        records.append({'file': path.name, 'bytes': path.stat().st_size,
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    (directory/'SHA256SUMS').write_text(''.join(f"{r['sha256']}  {r['file']}\n" for r in records))
    manifest = {'version': version, 'verified': True, 'artifacts': records,
                'revision': os.environ.get('GITHUB_SHA', 'local build'),
                'module_count': len(source_modules), 'source_files_checked': len(required)}
    (directory/'release-manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('dist'))
    args = parser.parse_args()
    print(json.dumps(check_distributions(Path(__file__).resolve().parents[1], args.directory), sort_keys=True))


if __name__ == '__main__':
    main()
