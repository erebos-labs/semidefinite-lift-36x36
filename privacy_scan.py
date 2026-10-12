"""Deterministic public-file privacy gate; manual review remains a separate evidence level."""
if not __debug__:
    raise RuntimeError('Run without -O: exact assertions must remain enabled')

import hashlib,json,pathlib,re,sys
ROOT=pathlib.Path(__file__).resolve().parent
PATTERNS={
    'private-filesystem-root':r'/(?:Users|home|private|tmp|var/folders|opt/homebrew|Library/Frameworks)/[A-Za-z0-9_.-][^\s"\x27,;)\]}]*',
    'personal-library-identifier':r'(?:libfile|libdir)_[0-9a-f]{20,}',
    'personal-upload-identifier':r'file_[0-9a-f]{20,}',
    'private-conversation-link':r'(?:chatgpt[.]com|chat[.]openai[.]com)/c/[0-9a-f-]+',
    'private-library-link':r'oai-library' + r'://',
    'email-address':r'\b[A-Za-z0-9_.+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b',
    'private-key-material':r'BEGIN (?:[A-Z]+ )?PRIVATE KEY',
    'credential-format':r'\b(?:gh[pousr]_[A-Za-z0-9]{20,}|sk-(?:proj-)?[A-Za-z0-9_-]{20,})',
    'signed-storage-parameters':r'[?&](?:sig|X-Amz-Signature|X-Goog-Signature)=',
}
IGNORE={'__pycache__','.venv','.git'}
def scan():
    hits=[];count=0;files=[]
    for p in sorted(ROOT.rglob('*')):
        if IGNORE.intersection(p.relative_to(ROOT).parts):continue
        if p.is_symlink():raise RuntimeError('Public candidate contains a symlink')
        if not p.is_file():continue
        if p.suffix in {'.zip','.sqlite','.db','.png','.jpg','.jpeg','.gif','.pdf'}:
            raise RuntimeError('Unreviewed binary file in public candidate: '+str(p.relative_to(ROOT)))
        data=p.read_bytes()
        try:text=data.decode('utf-8')
        except UnicodeDecodeError:raise RuntimeError('Non-text file in public candidate: '+str(p.relative_to(ROOT)))
        for label,pattern in PATTERNS.items():
            for found in re.finditer(pattern,text):
                hits.append({'file':str(p.relative_to(ROOT)),'rule':label,'line':text.count('\n',0,found.start())+1})
        count+=1;files.append({'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(data).hexdigest()})
    return {'status':'FAIL' if hits else 'PASS','files_scanned':count,'hits':hits,'scope':'UTF-8 public files; no private directory inventory, account data, or network calls.','limits':'Pattern scan plus a separate manual content review; no claim that heuristic scanning proves absence of every possible secret.','files':files}
if __name__=='__main__':
    result=scan();print(json.dumps({k:v for k,v in result.items() if k!='files'},indent=2))
    raise SystemExit(1 if result['hits'] else 0)
