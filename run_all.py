#!/usr/bin/env python3
"""Run the complete exact v1.0.0 verification and save attributable logs.

Requires Python >=3.10 and SymPy (see requirements.txt). Run without -O.
The coefficient generator and verify.py alone use only the standard library.
"""
import datetime
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise RuntimeError("Run without -O: exact assertions must remain enabled")
    import sympy
    logs=ROOT/"tests"
    logs.mkdir(exist_ok=True)
    stages=[
        ("core",["verify.py"]),
        ("independent",["audits/independent_coefficients.py",
                        "--primary-coefficients","coefficients.json"]),
        ("matching",["audits/check_matching.py"]),
        ("determinant",["audits/check_determinant_identity.py"]),
    ]
    result={"version":"1.0.0","started_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "python":platform.python_version(),"sympy":sympy.__version__,
            "arithmetic":"exact rational/quadratic-field/polynomial identities",
            "stages":[],"status":"RUNNING"}
    for name,arguments in stages:
        proc=subprocess.run([sys.executable,*arguments],cwd=ROOT,text=True,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        log_path=logs/(name+".log")
        log_path.write_text(proc.stdout)
        entry={"name":name,"command":["python",*arguments],"exit_code":proc.returncode,
               "log":str(log_path.relative_to(ROOT)),
               "log_sha256":hashlib.sha256(log_path.read_bytes()).hexdigest()}
        result["stages"].append(entry)
        print(("PASS " if proc.returncode==0 else "FAIL ")+name,flush=True)
        if proc.returncode:
            result["status"]="FAIL"
            (logs/"run_summary.json").write_text(json.dumps(result,indent=2)+"\n")
            print(proc.stdout)
            return proc.returncode
    result["status"]="PASS"
    result["completed_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
    result["scope_limits"]=[
        "Written proof supplies universal cone equivalence and logical deductions",
        "No Lean formalization or independent human peer review is claimed",
        "No numerical conditioning/performance or global minimality claim",
        "The independent stage repeats verification of the same fixed audit, not another search"]
    (logs/"run_summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS all exact verification stages; see tests/run_summary.json")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
