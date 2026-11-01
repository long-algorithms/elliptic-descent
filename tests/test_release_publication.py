"""Draft publication protocol tests; not substitutes for mathematical audits."""
import copy
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
import publish_release as publisher
sys.path.pop(0)
SHA = 'a'*40


def draft_release(sha=SHA):
    return dict(id=123, tag_name='v0.5.0', target_commitish=sha, draft=True,
                prerelease=False, published_at=None, assets=[])


def test_draft_lookup_uses_authenticated_listing_not_tag_endpoint(monkeypatch):
    draft = draft_release()
    seen = []
    def fake_api(path, **kwargs):
        seen.append(path)
        assert '/releases?' in path and '/tags/' not in path
        return [{'id': i+1000, 'tag_name': f'other-{i}'} for i in range(100)] if path.endswith('&page=1') else [draft]
    monkeypatch.setattr(publisher, 'api', fake_api)
    assert publisher.find_release('owner/repo', 'v0.5.0') == draft
    assert len(seen) == 2
    monkeypatch.setattr(publisher, 'api', lambda *a, **k: [])
    assert publisher.find_release('owner/repo', 'v0.5.0') is None


def test_duplicate_draft_tags_across_pages_are_rejected(monkeypatch):
    draft = draft_release()
    def fake_api(path, **kwargs):
        return [draft]+[dict(tag_name=f'other-{i}') for i in range(99)] if path.endswith('&page=1') else [{**draft, 'id':456}]
    monkeypatch.setattr(publisher, 'api', fake_api)
    with pytest.raises(ValueError, match='multiple releases'):
        publisher.find_release('owner/repo', 'v0.5.0')


def test_api_patch_uses_json_stdin_preserving_boolean_types(monkeypatch):
    import json
    def fake_gh(*args, **kwargs):
        assert args == ('api','repos/owner/repo/releases/123','--method','PATCH','--input','-')
        assert json.loads(kwargs['input_data']) == {'draft':False, 'make_latest':'true'}
        return b'{"id":123}'
    monkeypatch.setattr(publisher, 'gh', fake_gh)
    assert publisher.api('repos/owner/repo/releases/123',method='PATCH',
                         payload={'draft':False,'make_latest':'true'}) == {'id':123}


def test_explicit_draft_recovery_requires_an_untagged_ancestor(monkeypatch):
    old = 'b'*40
    draft = draft_release(old)
    monkeypatch.setattr(publisher, 'tag_commit', lambda *a: None)
    monkeypatch.setattr(publisher, 'api', lambda *a, **k: dict(status='ahead', merge_base_commit={'sha':old}))
    publisher.replacement_allowed('owner/repo',draft,'v0.5.0',SHA,'123@'+old)
    monkeypatch.setattr(publisher, 'tag_commit', lambda *a: old)
    with pytest.raises(ValueError, match='existing release tag'):
        publisher.replacement_allowed('owner/repo',draft,'v0.5.0',SHA,'123@'+old)
    monkeypatch.setattr(publisher, 'tag_commit', lambda *a: None)
    monkeypatch.setattr(publisher, 'api', lambda *a, **k: dict(status='diverged', merge_base_commit={'sha':'c'*40}))
    with pytest.raises(ValueError, match='descendant'):
        publisher.replacement_allowed('owner/repo',draft,'v0.5.0',SHA,'123@'+old)


@pytest.mark.parametrize('replacement', [None, '', '123', '456@'+'b'*40, '123@'+'c'*40])
def test_unknown_draft_cannot_be_retargeted(monkeypatch,replacement):
    monkeypatch.setattr(publisher, 'tag_commit', lambda *a: None)
    with pytest.raises(ValueError):
        publisher.replacement_allowed('owner/repo',draft_release('b'*40),'v0.5.0',SHA,replacement)


@pytest.mark.parametrize('field,value', [('draft',False),('published_at','2026-11-01T00:00:00Z'),
    ('prerelease',True),('tag_name','v0.6.0'),('id',False)])
def test_published_or_altered_draft_is_never_recovered(monkeypatch,field,value):
    draft = {**draft_release('b'*40),field:value}
    with pytest.raises(ValueError):
        publisher.replacement_allowed('owner/repo',draft,'v0.5.0',SHA,'123@'+'b'*40)


def simulated_release(monkeypatch,tmp_path,*,old=False,empty=False,corrupt=None,
                      move_main=False,move_tag=False,placeholder=False):
    """Publication metadata only. Omitted tag_name reproduces the real failure."""
    import hashlib
    files = tmp_path/'output'
    files.mkdir()
    (files/'example.whl').write_bytes(b'audited distribution')
    (files/'SHA256SUMS').write_bytes(b'audited checksum file')
    state = dict(release=None if empty else draft_release('b'*40 if old else SHA), calls=[],
                 tag=None, main_reads=0, tag_reads=0, next_asset=1000)
    if placeholder:
        state['release']['tag_name'] = 'untagged-observed'
    def fake_api(path, optional=False, *, method=None, payload=None):
        state['calls'].append((path,method,payload))
        if path.endswith('/git/ref/heads/main'):
            state['main_reads'] += 1
            return {'object':{'sha':'c'*40 if move_main and state['main_reads']>1 else SHA}}
        if '/compare/' in path:
            return {'status':'ahead','merge_base_commit':{'sha':'b'*40}}
        if '/releases?' in path:
            return [] if state['release'] is None else [copy.deepcopy(state['release'])]
        if path == 'repos/owner/repo/releases' and method == 'POST':
            assert state['release'] is None
            state['release'] = {**draft_release(), **payload}
            return copy.deepcopy(state['release'])
        if path.startswith('repos/owner/repo/releases/123/assets?'):
            return copy.deepcopy(state['release']['assets'])
        # Neither release-by-tag nor a different numeric identity may be used.
        assert path == 'repos/owner/repo/releases/123', path
        if method == 'PATCH':
            assert state['release']['draft'] and state['tag'] is None
            # PR18's omitted tag_name really produced an untagged placeholder.
            state['release']['tag_name'] = payload.get('tag_name', 'untagged-observed')
            state['release'].update(payload)
            if payload.get('draft') is False:
                state['release']['published_at'] = '2026-11-01T00:00:00Z'
                state['tag'] = payload['target_commitish']
        return copy.deepcopy(state['release'])
    def fake_gh(*args, **kwargs):
        from urllib.parse import parse_qs, urlparse
        state['calls'].append(args)
        assert args[0] == 'api'  # No CLI create/upload by tag for a draft.
        if args[2:4] == ('--method', 'DELETE'):
            identifier = int(args[1].rsplit('/',1)[1])
            state['release']['assets'] = [a for a in state['release']['assets'] if a['id'] != identifier]
            return b''
        assert args[1].startswith('https://uploads.github.com/repos/owner/repo/releases/123/assets?')
        assert args[2:6] == ('--method','POST','--header','Content-Type: application/octet-stream')
        assert args[6] == '--input'
        assert state['release']['target_commitish'] == SHA and state['release']['draft']
        p = Path(args[7])
        assert parse_qs(urlparse(args[1]).query)['name'] == [p.name]
        state['next_asset'] += 1
        state['release']['assets'].append(dict(id=state['next_asset'],name=p.name,state='uploaded',
            size=p.stat().st_size,digest='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()))
        if corrupt == 'checksum':
            state['release']['assets'][0]['digest'] = 'sha256:'+'0'*64
        elif corrupt == 'duplicate':
            state['release']['assets'].append(state['release']['assets'][0])
        elif corrupt == 'missing':
            state['release']['assets'].pop()
        elif corrupt == 'state':
            state['release']['prerelease'] = True
        elif corrupt == 'identity':
            state['release']['tag_name'] = 'unrelated'
        return b'uploaded'
    def fake_tag(*args):
        state['tag_reads'] += 1
        return 'c'*40 if move_tag and state['tag_reads']>1 else state['tag']
    monkeypatch.setattr(publisher,'api',fake_api)
    monkeypatch.setattr(publisher,'gh',fake_gh)
    monkeypatch.setattr(publisher,'tag_commit',fake_tag)
    return files,state


@pytest.mark.parametrize('mode',['new','retry','recover'])
def test_publish_lifecycle_works_when_tag_endpoint_cannot_see_drafts(monkeypatch,tmp_path,mode):
    files,state = simulated_release(monkeypatch,tmp_path,old=mode=='recover',empty=mode=='new')
    initial = copy.deepcopy(state['release'])
    final = publisher.stage_and_publish('owner/repo',SHA,'v0.5.0',files,initial,
                                        '123@'+'b'*40 if mode=='recover' else None)
    assert final['draft'] is False and final['target_commitish'] == SHA
    assert state['tag'] == SHA
    writes = [c for c in state['calls'] if len(c)==3 and c[1]=='PATCH']
    assert len(writes) == (2 if mode=='recover' else 1)
    assert all(c[2]['tag_name'] == 'v0.5.0' and c[2]['target_commitish'] == SHA for c in writes)
    assert writes[-1][2] == dict(tag_name='v0.5.0', target_commitish=SHA, draft=False,
                                  prerelease=False, make_latest='true')


@pytest.mark.parametrize('corrupt',['checksum','duplicate','missing','state','identity'])
def test_bad_uploaded_assets_or_changed_release_block_publication(monkeypatch,tmp_path,corrupt):
    files,state = simulated_release(monkeypatch,tmp_path,corrupt=corrupt)
    with pytest.raises(ValueError):
        publisher.stage_and_publish('owner/repo',SHA,'v0.5.0',files,copy.deepcopy(state['release']))
    assert not any(len(c)==3 and c[1]=='PATCH' for c in state['calls'])
    assert state['tag'] is None


@pytest.mark.parametrize('mutation',['main','tag'])
def test_moving_ref_during_upload_blocks_publication(monkeypatch,tmp_path,mutation):
    files,state = simulated_release(monkeypatch,tmp_path,move_main=mutation=='main',move_tag=mutation=='tag')
    with pytest.raises(ValueError):
        publisher.stage_and_publish('owner/repo',SHA,'v0.5.0',files,copy.deepcopy(state['release']))
    assert not any(len(c)==3 and c[1]=='PATCH' for c in state['calls'])


def test_failed_mathematical_inventory_prevents_any_draft_mutation(monkeypatch):
    monkeypatch.setenv('GITHUB_REPOSITORY','owner/repo')
    monkeypatch.setattr(publisher.subprocess,'check_output',lambda *a,**k: SHA+'\n')
    old = 'b'*40
    def fake_api(path, **kwargs):
        if '/git/ref/heads/main' in path:
            return {'object':{'sha':SHA}}
        if '/compare/' in path:
            return {'status':'ahead','merge_base_commit':{'sha':old}}
        pytest.fail('unexpected API mutation or request: '+path)
    monkeypatch.setattr(publisher,'api',fake_api)
    monkeypatch.setattr(publisher,'find_release',lambda *a: draft_release(old))
    monkeypatch.setattr(publisher,'tag_commit',lambda *a: None)
    monkeypatch.setattr(publisher,'pages',lambda *a: [])
    monkeypatch.setattr(publisher,'selected_runs',lambda *a: {'test':{'id':1}})
    def bad_inventory(*args):
        raise ValueError('missing genuine mathematical evidence')
    monkeypatch.setattr(publisher,'audit_directory',bad_inventory)
    monkeypatch.setattr(publisher,'stage_and_publish',lambda *a,**k: pytest.fail('mutation before audit'))
    with pytest.raises(ValueError, match='mathematical evidence'):
        publisher.publish('owner/repo',SHA,replacement='123@'+old)


def test_every_release_gate_runs_on_main_without_path_filtering():
    # A tooling/documentation-only merge must not strand the exact-head gate.
    import re
    from release_gate import WORKFLOWS
    for name in WORKFLOWS:
        text = (ROOT/'.github/workflows'/name).read_text()
        match = re.search(r'^  push:\n((?:    .*\n)*)', text, re.M)
        assert match is not None, name
        block = match[1]
        assert 'paths:' not in block and 'paths-ignore:' not in block, name
        assert 'branches-ignore:' not in block, name
        assert block.strip() == 'branches: [main]', name


def test_group_replay_has_separate_budget_without_omitting_any_fixture():
    import re
    from release_gate import GROUPS
    text = (ROOT/'.github/workflows/group-benchmarks.yml').read_text()
    assert '  replay:\n' in text and '    needs: group\n' in text
    assert text.count('    timeout-minutes: 20\n') == 2
    matrices = re.findall(r'fixture: \[(.*?)\]', text)
    assert len(matrices) == 2
    assert all(set(v.strip() for v in m.split(',')) == set(GROUPS) for m in matrices)
    assert 'benchmarks/four_selmer_group.py --fixture ${{ matrix.fixture }} --verify' in text
    assert 'benchmarks/replay_four_group.py four-selmer-group-${{ matrix.fixture }}.jsonl' in text
    assert 'name: group-construction-${{ matrix.fixture }}' in text
    assert 'construction-${{ matrix.fixture }}.record' in text
    assert 'name: group-${{ matrix.fixture }}' in text
    assert 'uses: actions/download-artifact@v4' in text


def test_observed_placeholder_recovered_by_explicit_numeric_identity(monkeypatch,tmp_path):
    files,state = simulated_release(monkeypatch,tmp_path,old=True,placeholder=True)
    old = 'b'*40
    original = copy.deepcopy(state['release'])
    found = publisher.locate_release('owner/repo','v0.5.0',SHA,'123@'+old,'untagged-observed')
    assert found == original
    final = publisher.stage_and_publish('owner/repo',SHA,'v0.5.0',files,found,
                                         '123@'+old,'untagged-observed')
    assert final['id'] == 123 and final['tag_name'] == 'v0.5.0'
    assert final['draft'] is False and state['tag'] == SHA
    assert not any(c[0] == 'repos/owner/repo/releases' and c[1] == 'POST'
                   for c in state['calls'] if len(c) == 3)


@pytest.mark.parametrize('change', [{'id':456}, {'target_commitish':'c'*40},
    {'tag_name':'untagged-another'}, {'draft':False}, {'prerelease':True},
    {'published_at':'2026-11-01T00:00:00Z'}])
def test_placeholder_recovery_rejects_any_identity_or_state_mismatch(monkeypatch,change):
    draft = {**draft_release('b'*40), 'tag_name':'untagged-observed', **change}
    monkeypatch.setattr(publisher,'find_release',lambda *a:None)
    monkeypatch.setattr(publisher,'api',lambda *a,**k:draft)
    monkeypatch.setattr(publisher,'tag_commit',lambda *a:None)
    with pytest.raises(ValueError):
        publisher.locate_release('owner/repo','v0.5.0',SHA,'123@'+'b'*40,'untagged-observed')


def test_placeholder_name_alone_never_authorizes_recovery(monkeypatch):
    monkeypatch.setattr(publisher,'find_release',lambda *a:None)
    monkeypatch.setattr(publisher,'api',lambda *a,**k:pytest.fail('request without explicit ID'))
    with pytest.raises(ValueError):
        publisher.locate_release('owner/repo','v0.5.0',SHA,recovery_tag='untagged-observed')


def test_asset_upload_by_id_reuses_identical_files_and_replaces_changed_expected_files(monkeypatch,tmp_path):
    import hashlib
    files,state = simulated_release(monkeypatch,tmp_path)
    original = dict(id=300,name='example.whl',state='uploaded',size=1,digest='sha256:'+'0'*64)
    same = dict(id=301,name='SHA256SUMS',state='uploaded',size=(files/'SHA256SUMS').stat().st_size,
                digest='sha256:'+hashlib.sha256((files/'SHA256SUMS').read_bytes()).hexdigest())
    state['release']['assets'] = [original,same]
    publisher.upload_assets_by_id('owner/repo',123,files)
    publisher.check_uploaded_assets(state['release'],files)
    assert any(c[:4] == ('api','repos/owner/repo/releases/assets/300','--method','DELETE') for c in state['calls'])
    assert not any(c[:4] == ('api','repos/owner/repo/releases/assets/301','--method','DELETE') for c in state['calls'])
    state['calls'].clear()
    publisher.upload_assets_by_id('owner/repo',123,files)
    assert all(len(c) == 3 and c[1] is None for c in state['calls'])


@pytest.mark.parametrize('names',[['unknown'],['example.whl','example.whl']])
def test_unknown_or_duplicate_assets_are_not_deleted(monkeypatch,tmp_path,names):
    files,state = simulated_release(monkeypatch,tmp_path)
    state['release']['assets'] = [dict(id=i+1,name=n) for i,n in enumerate(names)]
    with pytest.raises(ValueError):
        publisher.upload_assets_by_id('owner/repo',123,files)
    assert not any(c[0] == 'api' for c in state['calls'])


def test_asset_listing_is_exhaustive(monkeypatch):
    monkeypatch.setattr(publisher,'api',lambda path: [{'id':i} for i in range(100)]
                        if path.endswith('page=1') else [{'id':101}])
    assert len(publisher.release_assets('owner/repo',123)) == 101
