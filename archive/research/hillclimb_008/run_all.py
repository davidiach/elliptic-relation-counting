"""Retain the exact readout audit, certificates, and source hashes."""
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
    paths=list(ROOT.glob('*.py'))+list(ROOT.glob('*.tex'))
    paths.extend([ROOT.parent/'hillclimb_007/cluster_audit.py',
                  ROOT.parent/'hillclimb_007/cluster_bounds.tex',
                  ROOT.parent/'hillclimb_006/nonuniform_grid.py'])
    return {str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    run=ROOT/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    run.mkdir(parents=True)
    before=hashes();start=time.monotonic()
    env=dict(os.environ,PYTHONHASHSEED='0',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONOPTIMIZE',None)
    stdout_path=run/'stdout.txt';stderr_path=run/'stderr.txt'
    with stdout_path.open('w',encoding='utf-8') as out,stderr_path.open('w',encoding='utf-8') as err:
        try:
            process=subprocess.run([sys.executable,'-B',str(ROOT/'readout_audit.py')],
                                   cwd=ROOT,env=env,text=True,stdout=out,stderr=err,timeout=240)
            code=process.returncode
        except subprocess.TimeoutExpired:
            code=124
            err.write('Audit exceeded its explicit 240-second process limit.\n')
    summary=dict(returncode=code,seconds=round(time.monotonic()-start,3),
                 source_sha256=before,sources_unchanged=before==hashes(),python=sys.version.split()[0])
    summary['status']='failed'
    if code==0:
        data=json.loads((ROOT/'readout_results.json').read_text(encoding='utf-8'))
        shutil.copyfile(ROOT/'readout_results.json',run/'readout_results.json')
        summary['result_sha256']=hashlib.sha256((run/'readout_results.json').read_bytes()).hexdigest()
        if summary['sources_unchanged'] and data['status']=='passed':summary['status']='passed'
        summary['counts']=data['counts']
    (run/'execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (ROOT/'latest_run.json').write_text(json.dumps(dict(path=str(run.relative_to(ROOT)),status=summary['status']),indent=2)+'\n',encoding='utf-8')
    print(stdout_path.read_text(encoding='utf-8'),end='')
    errors=stderr_path.read_text(encoding='utf-8')
    if errors:print(errors,file=sys.stderr,end='')
    print(json.dumps(dict(status=summary['status'],run=str(run))))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
