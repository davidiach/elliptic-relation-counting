"""Retain actual-packet rank and parity certificates with source identities."""
from datetime import datetime,timezone
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
    paths=list(ROOT.glob('*.py'))+list(ROOT.glob('*.tex'))+[ROOT/'basis_labels.json']
    for name in ('hillclimb_003/c38_slices.py','hillclimb_002/c38_network.py',
                 'hillclimb_001/c38_simplex.py','hillclimb_001/common.py',
                 'hillclimb_007/cluster_audit.py','hillclimb_006/nonuniform_grid.py',
                 'hillclimb_008/optimal_readout.py','hillclimb_008/verify_certificate.py',
                 'hillclimb_008/readout_results.json'):
        paths.append(ROOT.parent/name)
    paths.extend((ROOT.parent/'hillclimb_001/vendor').glob('*.py'))
    return {str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    run=ROOT/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    run.mkdir(parents=True);before=hashes();start=time.monotonic()
    env=dict(os.environ,PYTHONHASHSEED='0',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONOPTIMIZE',None)
    with (run/'stdout.txt').open('w',encoding='utf-8') as out,(run/'stderr.txt').open('w',encoding='utf-8') as err:
        try:
            code=subprocess.run([sys.executable,'-B',str(ROOT/'space_audit.py')],cwd=ROOT,
                                env=env,text=True,stdout=out,stderr=err,timeout=180).returncode
        except subprocess.TimeoutExpired:
            code=124;err.write('Audit exceeded its explicit 180-second process limit.\n')
    summary=dict(status='failed',returncode=code,seconds=round(time.monotonic()-start,3),
                 source_sha256=before,sources_unchanged=before==hashes(),python=sys.version.split()[0])
    if code==0:
        data=json.loads((ROOT/'space_results.json').read_text(encoding='utf-8'))
        shutil.copyfile(ROOT/'space_results.json',run/'space_results.json')
        summary['result_sha256']=hashlib.sha256((run/'space_results.json').read_bytes()).hexdigest()
        if summary['sources_unchanged'] and data['status']=='passed':summary['status']='passed'
    (run/'execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (ROOT/'latest_run.json').write_text(json.dumps(dict(path=str(run.relative_to(ROOT)),status=summary['status']),indent=2)+'\n',encoding='utf-8')
    print((run/'stdout.txt').read_text(encoding='utf-8'),end='')
    error=(run/'stderr.txt').read_text(encoding='utf-8')
    if error:print(error,file=sys.stderr,end='')
    print(json.dumps(dict(status=summary['status'],run=str(run))))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
