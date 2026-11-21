"""Source-hashed checkpoint for the compression and refactoring iteration."""
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

ROOT = Path(__file__).resolve().parent
JOBS = (('probe_compression.py','compression_results.json'),
        ('verify_factoring.py','factoring_results.json'))


def hashes():
    paths = list(ROOT.glob('*.py'))+list(ROOT.glob('*.tex'))
    for name in ('hillclimb_001','hillclimb_002'):
        previous = ROOT.parent/name
        paths += list(previous.glob('*.py'))
        paths += list((previous/'vendor').glob('*.py'))
    paths += [ROOT.parent/'hillclimb_001/source_manifest.json']
    return {str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    run = ROOT/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    run.mkdir(parents=True)
    before = hashes()
    env = dict(os.environ,PYTHONHASHSEED='0',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONOPTIMIZE',None)
    records = []
    for script,result in JOBS:
        start = time.monotonic()
        process = subprocess.run([sys.executable,'-B',str(ROOT/script)],cwd=ROOT,
                                 env=env,text=True,capture_output=True,timeout=300)
        (run/(script+'.stdout.txt')).write_text(process.stdout,encoding='utf-8')
        (run/(script+'.stderr.txt')).write_text(process.stderr,encoding='utf-8')
        record = dict(script=script,returncode=process.returncode,seconds=round(time.monotonic()-start,3),result=result)
        if process.returncode == 0:
            record['status'] = json.loads((ROOT/result).read_text())['status']
            shutil.copyfile(ROOT/result,run/result)
        else:
            record['status'] = 'failed'
        records.append(record)
        print(json.dumps(record),flush=True)
    unchanged = before == hashes()
    passed = unchanged and all(r['status'] == 'passed' for r in records)
    summary = dict(status='passed' if passed else 'failed',python=sys.version.split()[0],records=records,
                   source_hashes=before,sources_unchanged=unchanged,
                   scope='Finite actual C38 rank and projection probes, exact factored joins for a restricted family, and complete graph refactoring checks. No improved ECDLP runtime or full recursive compact-representation closure is certified.')
    (run/'execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (ROOT/'latest_run.json').write_text(json.dumps({'path':str(run.relative_to(ROOT)),'status':summary['status']},indent=2)+'\n')
    print(json.dumps({'status':summary['status'],'run':str(run)}))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
