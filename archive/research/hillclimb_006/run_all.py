"""Keep a source-hashed execution checkpoint for the nonuniform-grid audit."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

ROOT=Path(__file__).resolve().parent


def hashes():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(list(ROOT.glob('*.py'))+list(ROOT.glob('*.tex')))}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    run=ROOT/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    run.mkdir(parents=True)
    before=hashes()
    env=dict(os.environ,PYTHONHASHSEED='0',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONOPTIMIZE',None)
    start=time.monotonic()
    process=subprocess.run([sys.executable,'-B',str(ROOT/'nonuniform_grid.py')],cwd=ROOT,
                           env=env,text=True,capture_output=True,timeout=180)
    (run/'stdout.txt').write_text(process.stdout,encoding='utf-8')
    (run/'stderr.txt').write_text(process.stderr,encoding='utf-8')
    summary=dict(returncode=process.returncode,seconds=round(time.monotonic()-start,3),
                 source_sha256=before,sources_unchanged=before==hashes(),python=sys.version.split()[0])
    if process.returncode==0:
        data=json.loads((ROOT/'grid_results.json').read_text())
        shutil.copyfile(ROOT/'grid_results.json',run/'grid_results.json')
        summary['result_sha256']=hashlib.sha256((run/'grid_results.json').read_bytes()).hexdigest()
        summary['status']='passed' if summary['sources_unchanged'] and data['status']=='passed' else 'failed'
    else:
        summary['status']='failed'
    (run/'execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (ROOT/'latest_run.json').write_text(json.dumps(dict(path=str(run.relative_to(ROOT)),status=summary['status']),indent=2)+'\n')
    print(process.stdout,end='')
    if process.stderr:
        print(process.stderr,file=sys.stderr,end='')
    print(json.dumps(dict(status=summary['status'],run=str(run))))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':
    raise SystemExit(main())
