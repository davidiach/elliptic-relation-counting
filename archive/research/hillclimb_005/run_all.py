"""Retain source-hashed executions of the balanced-packet research probes."""
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
JOBS=(('balanced_projection.py','balanced_results.json'),
      ('normalized_permanent.py','normalized_results.json'))


def hashes():
    paths=list(ROOT.glob('*.py'))+list(ROOT.glob('*.tex'))
    for name in ('hillclimb_001','hillclimb_002','hillclimb_003'):
        folder=ROOT.parent/name
        paths+=list(folder.glob('*.py'))+list((folder/'vendor').glob('*.py'))
    paths+=[ROOT.parent/'hillclimb_001/source_manifest.json']
    return {str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    run=ROOT/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    run.mkdir(parents=True)
    before=hashes()
    env=dict(os.environ,PYTHONHASHSEED='0',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONOPTIMIZE',None)
    records=[]
    for script,record in JOBS:
        start=time.monotonic()
        process=subprocess.run([sys.executable,'-B',str(ROOT/script)],cwd=ROOT,
                               env=env,text=True,capture_output=True,timeout=180)
        (run/(script+'.stdout.txt')).write_text(process.stdout,encoding='utf-8')
        (run/(script+'.stderr.txt')).write_text(process.stderr,encoding='utf-8')
        item=dict(script=script,returncode=process.returncode,seconds=round(time.monotonic()-start,3))
        if process.returncode==0:
            item['status']=json.loads((ROOT/record).read_text())['status']
            shutil.copyfile(ROOT/record,run/record)
            item['result_sha256']=hashlib.sha256((run/record).read_bytes()).hexdigest()
        else:
            item['status']='failed'
        records.append(item)
        print(json.dumps(item),flush=True)
    unchanged=before==hashes()
    passed=unchanged and all(r['status']=='passed' for r in records)
    summary=dict(status='passed' if passed else 'failed',records=records,source_sha256=before,
                 sources_unchanged=unchanged,python=sys.version.split()[0],
                 scope='Finite balanced C38 projection identity and counterexamples to specific compression shortcuts. No complete ECDLP improvement.')
    (run/'execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (ROOT/'latest_run.json').write_text(json.dumps(dict(path=str(run.relative_to(ROOT)),status=summary['status']),indent=2)+'\n')
    print(json.dumps(dict(status=summary['status'],run=str(run))))
    return 0 if passed else 1


if __name__=='__main__':
    raise SystemExit(main())
