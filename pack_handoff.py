#!/usr/bin/env python3
"""Package this completed handoff with SHA-256 inventory and exact ZIP check.

Run only after verification and document review. This performs no network or
repository operations. Generated caches and virtual environments are excluded.
"""
if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

from pathlib import Path
import hashlib
import json
import sys
import subprocess
import zipfile
from release_utils import sanitize_full_replay_logs
from privacy_scan import scan


def main():
    root=Path(__file__).resolve().parent
    version=(root/"VERSION").read_text().strip()
    if version!="1.0.2":
        raise ValueError("Update the packer intentionally for another release")
    # Recompute all four stages; a saved PASS summary cannot authorize packaging.
    subprocess.run([sys.executable, str(root/"run_all.py")], cwd=root, check=True)
    summary=json.loads((root/"tests/run_summary.json").read_text())
    if summary.get("status")!="PASS" or len(summary.get("stages",[]))!=4:
        raise ValueError("A completed four-stage replay is required before packaging")
    sanitize_full_replay_logs()
    privacy=scan()
    if privacy['hits']:
        raise RuntimeError('Public privacy scan failed: '+json.dumps(privacy['hits']))
    ignored={"__pycache__",".venv",".git"}
    files=sorted(p for p in root.rglob("*") if p.is_file()
                 and not ignored.intersection(p.relative_to(root).parts)
                 and p.name!="SHA256SUMS.txt")
    manifest=root/"SHA256SUMS.txt"
    manifest.write_text("".join(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+
                                p.relative_to(root).as_posix()+"\n" for p in files))
    files=sorted(files+[manifest])
    archive=root.parent/f"semidefinite_lift_public_v{version}.zip"
    prefix=f"semidefinite_lift_public_v{version}"
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in files:
            info=zipfile.ZipInfo(prefix+"/"+p.relative_to(root).as_posix(),
                                 date_time=(2026,10,12,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o100644<<16
            z.writestr(info,p.read_bytes())
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ValueError("ZIP CRC failure")
        for p in files:
            if z.read(prefix+"/"+p.relative_to(root).as_posix())!=p.read_bytes():
                raise ValueError("ZIP byte mismatch: "+str(p))
    print(json.dumps({"archive":archive.name,"version":version,"files":len(files),
                      "bytes":archive.stat().st_size,
                      "sha256":hashlib.sha256(archive.read_bytes()).hexdigest(),
                      "zip_crc_and_all_entry_bytes":"PASS"},indent=2))


if __name__=="__main__":
    main()
