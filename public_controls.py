"""Scoped checker hardening acceptance tests. No new mathematical search."""
if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

import datetime, json, os, pathlib, re, shutil, subprocess, sys, tempfile
ROOT=pathlib.Path(__file__).resolve().parent
LOGS=ROOT/'evidence/sanitized/hardening';LOGS.mkdir(parents=True,exist_ok=True)
from release_utils import write_sanitized_log, sanitize_text
records=[]
def require(condition,message):
    if not condition:raise RuntimeError(message)
def run(label,arguments,cwd=ROOT,mode=None,success=True,guard=False):
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);flags=[]
    if mode in ('-O','-OO'):flags=[mode]
    elif mode:env['PYTHONOPTIMIZE']=mode.split('=')[1]
    p=subprocess.run([sys.executable,*flags,*arguments],cwd=cwd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
    log=LOGS/(label+'.log');log_metadata=write_sanitized_log(log,p.stdout,label)
    require((p.returncode==0)==success, label+' unexpected exit '+str(p.returncode))
    if guard:
        require('exact assertions must remain enabled' in p.stdout,label+' did not hit optimization guard')
        require(re.search(r'^\s*PASS\b|"status"\s*:\s*"PASS"',p.stdout,re.M) is None,label+' emitted a PASS claim')
    records.append({'control':label,'command':arguments,'mode':mode or 'normal','expected_success':success,'exit_code':p.returncode,'expectation_met':True,'no_PASS_with_assertions_disabled':guard,'log':str(log.relative_to(ROOT)),'log_representation':log_metadata})
    print(label+' expected outcome reproduced',flush=True)
    return p
direct=[('core',['verify.py']),('independent',['audits/independent_coefficients.py','--primary-coefficients','coefficients.json']),('matching',['audits/check_matching.py']),('determinant',['audits/check_determinant_identity.py']),('full',['run_all.py'])]
supplemental=[('third-coefficients',['audits/additional/audit_independent.py']),('direct-generic-matching',['audits/additional/matching_direct_audit.py']),('singular-roots',['audits/additional/singular_root_audit.py'])]
for name,args in direct+supplemental:run('normal-'+name,args)
guarded=direct+supplemental+[('packer',['pack_handoff.py']),('hardening-controls',['public_controls.py']),('privacy-scan',['privacy_scan.py'])]
for name,args in guarded:
    for mode in ['-O','-OO','PYTHONOPTIMIZE=1','PYTHONOPTIMIZE=2']:
        label='disabled-'+name+'-'+mode.replace('=','-').lstrip('-')
        run(label,args,mode=mode,success=False,guard=True)

with tempfile.TemporaryDirectory(prefix='semidefinite-v1.0.1-negative-') as temporary:
    temporary=pathlib.Path(temporary)
    def copy(label):
        d=temporary/label;shutil.copytree(ROOT,d,ignore=shutil.ignore_patterns('__pycache__','.venv','hardening-logs'));return d
    def corrupt_coeff(d):
        p=d/'coefficients.json';value=json.loads(p.read_text());value['coefficients'][0][0][2][0]+=1;p.write_text(json.dumps(value))
    d=copy('coefficient-corrupt');corrupt_coeff(d)
    run('corrupt-coefficient-core-rejected',['verify.py'],cwd=d,success=False)
    run('corrupt-coefficient-independent-rejected',['audits/independent_coefficients.py','--primary-coefficients','coefficients.json'],cwd=d,success=False)
    d=copy('missing-matrix');p=d/'coefficients.json';v=json.loads(p.read_text());v['coefficients'].pop();p.write_text(json.dumps(v))
    run('missing-matrix-rejected',['audits/independent_coefficients.py','--primary-coefficients','coefficients.json'],cwd=d,success=False)
    d=copy('matching-corrupt');p=d/'audits/check_matching.py';text=p.read_text();before='RF * RF.T - Q(y))';require(text.count(before)==1,'matching mutation anchor');p.write_text(text.replace(before,'RF * RF.T - Q(y) + s.eye(3))'))
    run('false-matching-normal-rejected',['audits/check_matching.py'],cwd=d,success=False)
    for mode in ['-O','-OO','PYTHONOPTIMIZE=1','PYTHONOPTIMIZE=2']:
        run('false-matching-disabled-'+mode.replace('=','-').lstrip('-'),['audits/check_matching.py'],cwd=d,mode=mode,success=False,guard=True)
    d=copy('determinant-corrupt');p=d/'audits/check_determinant_identity.py';text=p.read_text();before='check_zero(b.T*hat_adj*b-phi)';require(text.count(before)==1,'determinant mutation anchor');p.write_text(text.replace(before,'check_zero(b.T*hat_adj*b+phi)'))
    run('false-Phi-sign-rejected',['audits/check_determinant_identity.py'],cwd=d,success=False)
    d=copy('source-corrupt');p=d/'source/build/sections/angular.tex';p.write_text(p.read_text()+'\n% hardening negative control\n')
    run('source-corruption-rejected',['verify.py'],cwd=d,success=False)
    d=copy('forged-rank-certificate')
    for path in ['certificates/rank_kernel.json','certificates/verification.json','audits/redundancy_certificate.json']:
        (d/path).write_text(json.dumps({'status':'PASS','auxiliary_map_rank':0,'full_map_rank':0,'common_kernel_dimension':36}))
    run('forged-rank-certificate-recomputed',['verify.py'],cwd=d)
    rank=json.loads((d/'certificates/rank_kernel.json').read_text());v=json.loads((d/'certificates/verification.json').read_text())
    require(rank['full_map_rank']==47 and rank['auxiliary_map_rank']==24 and rank['common_kernel_dimension']==0,'saved false certificate trusted')
    require(v['coefficient_verification_complete'] is True and v['universal_matching_checker_included_here'] is False,'verification scope changed or obscured')
    (LOGS/'forged-rank-certificate-recomputed-output.json').write_text(json.dumps(rank,indent=2)+'\n')
    d=copy('forged-comparison');(d/'audits/independent_comparison.json').write_text(json.dumps({'status':'PASS','coefficient_matrices_compared':0}))
    run('forged-comparison-recomputed',['audits/independent_coefficients.py','--primary-coefficients','coefficients.json'],cwd=d)
    comparison=json.loads((d/'audits/independent_comparison.json').read_text())
    require(comparison['coefficient_matrices_compared']==47 and comparison['full_matrix_scalar_entries_compared']==60912,'saved comparison trusted')
    (LOGS/'forged-comparison-recomputed-output.json').write_text(json.dumps(comparison,indent=2)+'\n')
    d=copy('forged-PASS-summary');corrupt_coeff(d)
    (d/'tests/run_summary.json').write_text(json.dumps({'status':'PASS','stages':[{'exit_code':0}]*4}))
    run('forged-PASS-summary-cannot-authorize-replay',['run_all.py'],cwd=d,success=False)
    summary=json.loads((d/'tests/run_summary.json').read_text())
    require(summary['status']=='FAIL' and summary['stages'][0]['exit_code']!=0,'forged summary preserved')
    (LOGS/'forged-PASS-summary-replaced.json').write_text(json.dumps(summary,indent=2)+'\n')
    (d/'tests/run_summary.json').write_text(json.dumps({'status':'PASS','stages':[{'exit_code':0}]*4}))
    run('forged-PASS-summary-cannot-authorize-packaging',['pack_handoff.py'],cwd=d,success=False)
    # Synthetic leak controls: values are assembled at runtime so the detector's
    # own public source never contains a private-looking value.
    fixtures=[('filesystem','/'+'Users'+'/'+'example-private'+'/'+'fixture.txt'),
              ('library','libfile_'+'b'*32),
              ('email','person'+'@'+'example.invalid'),
              ('credential','ghp_'+'a'*40)]
    for label,content in fixtures:
        d=copy('privacy-'+label);(d/'deliberate-privacy-fixture.txt').write_text(content)
        p=run('privacy-fixture-'+label+'-rejected',['privacy_scan.py'],cwd=d,success=False)
        require(json.loads(p.stdout)['hits'],'Privacy scan was a vacuous oracle')
result={'status':'PASS','distribution_version':'1.0.2','mathematical_payload_version':'1.0.0','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sanitize_text(sys.version),'control_count':len(records),'normal_commands':len(direct+supplemental),'disabled_assertion_entrypoint_controls':4*len(guarded),'mutated_matching_disabled_controls':4,'records':records,'scope_limits':['The guards reject optimized execution; they do not implement assertion-free mathematical verification.','All original mathematical checks remain present.','Saved certificate outputs are recomputed, not treated as proof inputs.','Finite tests supplement the unchanged written proof; no new theorem review or research was performed.']}
(ROOT/'evidence/sanitized/hardening-summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
