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


def simulated_release(monkeypatch,tmp_path,*,old=False,empty=False,corrupt=None,move_main=False,move_tag=False):
    """Publication-protocol stub only; it supplies no mathematical evidence."""
    import hashlib
    files = tmp_path/'output'
    files.mkdir()
    (files/'example.whl').write_bytes(b'audited distribution')
    (files/'SHA256SUMS').write_bytes(b'audited checksum file')
    state = dict(release=None if empty else draft_release('b'*40 if old else SHA), calls=[],
                 tag=None, main_reads=0, tag_reads=0)
    def fake_api(path, optional=False, *, method=None, payload=None):
        state['calls'].append((path,method,payload))
        if path.endswith('/git/ref/heads/main'):
            state['main_reads'] += 1
            return {'object':{'sha':'c'*40 if move_main and state['main_reads']>1 else SHA}}
        if '/compare/' in path:
            return {'status':'ahead','merge_base_commit':{'sha':'b'*40}}
        if '/releases?' in path:
            return [] if state['release'] is None else [copy.deepcopy(state['release'])]
        # A release-by-tag request reproduces the real GitHub failure.
        assert path == 'repos/owner/repo/releases/123', path
        if method == 'PATCH':
            if 'target_commitish' in payload:
                assert state['release']['draft'] and state['tag'] is None
                state['release']['target_commitish'] = payload['target_commitish']
            if 'draft' in payload:
                assert payload == {'draft':False,'make_latest':'true'}
                state['release']['draft'] = False
                state['release']['published_at'] = '2026-11-01T00:00:00Z'
                state['tag'] = SHA
        return copy.deepcopy(state['release'])
    def fake_gh(*args, **kwargs):
        state['calls'].append(args)
        if args[:2] == ('release','create'):
            assert state['release'] is None
            state['release'] = draft_release()
            return b'created'
        assert args[:2] == ('release','upload')
        assert state['release']['target_commitish'] == SHA and state['release']['draft']
        state['release']['assets'] = [dict(name=p.name,state='uploaded',size=p.stat().st_size,
            digest='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()) for p in files.iterdir()]
        if corrupt == 'checksum':
            state['release']['assets'][0]['digest'] = 'sha256:'+'0'*64
        elif corrupt == 'duplicate':
            state['release']['assets'].append(state['release']['assets'][0])
        elif corrupt == 'missing':
            state['release']['assets'].pop()
        elif corrupt == 'state':
            state['release']['draft'] = False
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
    assert writes[-1][2] == {'draft':False,'make_latest':'true'}


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
