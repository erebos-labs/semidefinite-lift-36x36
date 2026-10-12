"""Release-log sanitization. Original execution logs are never called sanitized originals."""
import hashlib,json,os,pathlib,re

ROOT=pathlib.Path(__file__).resolve().parent
ABSOLUTE_PATH=re.compile(r'/(?:Users|home|private|tmp|var/folders|opt/homebrew|Library/Frameworks)/[^\s"\x27,;)\]}]+')

def sanitize_text(text):
    return ABSOLUTE_PATH.sub('<filesystem-path>',text)

def private_raw_dir():
    value=os.environ.get('SEMIDEFINITE_PRIVATE_RAW_DIR')
    if not value:return None
    path=pathlib.Path(value).resolve()
    if path==ROOT or ROOT in path.parents:
        raise RuntimeError('Raw log destination must be outside the public candidate')
    path.mkdir(parents=True,exist_ok=True)
    return path

def write_sanitized_log(destination,text,label):
    destination=pathlib.Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    raw=text.encode();digest=hashlib.sha256(raw).hexdigest();private=private_raw_dir()
    if private:
        (private/(label+'.raw.log')).write_bytes(raw)
    sanitized='SANITIZED EXECUTION LOG: filesystem paths redacted; outcomes and algebra unchanged.\n'+sanitize_text(text)
    destination.write_text(sanitized)
    return {'representation':'sanitized-copy','raw_sha256':digest,'sanitized_sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),'raw_retention':'external private directory' if private else 'not retained by this invocation'}

def sanitize_full_replay_logs():
    path=ROOT/'tests/run_summary.json';result=json.loads(path.read_text())
    for stage in result['stages']:
        log=ROOT/stage['log'];raw=log.read_text()
        meta=write_sanitized_log(log,raw,'full-replay-'+stage['name'])
        stage['original_runtime_log_sha256']=stage['log_sha256']
        stage['log_sha256']=meta['sanitized_sha256'];stage['log_representation']='sanitized-copy'
    result['log_representation']='sanitized-copy'
    result['sanitization_scope']='Filesystem roots only; raw outputs retained privately for preparation when configured.'
    path.write_text(json.dumps(result,indent=2)+'\n')
