"""Run this research iteration, preserving each run and its source hashes."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import uuid
import common


def main():
    root = Path(__file__).resolve().parent
    common.verify_sources()
    run = root/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    run.mkdir(parents=True)
    source_paths = sorted(root.glob('*.py'))
    source_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    env = dict(os.environ, PYTHONOPTIMIZE='0', PYTHONHASHSEED='0', PYTHONIOENCODING='utf-8',
               PYTHONDONTWRITEBYTECODE='1', HILLCLIMB_RESULTS=str(run))
    jobs = [('verify_c38_simplex.py','c38_simplex.json'),
            ('verify_adaptive_simplex.py','adaptive_simplex.json'),
            ('certify_bounds.py','bounds.json')]
    records = []
    for script, output in jobs:
        start = time.monotonic()
        proc = subprocess.run([sys.executable,'-B',script], cwd=root, env=env, capture_output=True, text=True, timeout=300)
        (run/(script+'.stdout.txt')).write_text(proc.stdout,encoding='utf-8')
        (run/(script+'.stderr.txt')).write_text(proc.stderr,encoding='utf-8')
        passed = proc.returncode == 0 and (run/output).exists()
        if passed:
            passed = json.loads((run/output).read_text())['status']=='passed'
        row = dict(script=script,returncode=proc.returncode,status='passed' if passed else 'failed',
                   seconds=round(time.monotonic()-start,3),result=output)
        records.append(row)
        print(json.dumps(row),flush=True)
        if not passed:
            print(proc.stderr)
    after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    unchanged = source_hashes == after
    common.verify_sources()
    status = 'passed' if unchanged and all(r['status']=='passed' for r in records) else 'failed'
    summary = dict(status=status,python=platform.python_version(),records=records,
                   source_hashes=source_hashes,sources_unchanged=unchanged,
                   scope='Finite checks and exact certificates supporting a conditional composition; not full proof certification.')
    (run/'execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (root/'latest_run.json').write_text(json.dumps({'path':str(run.relative_to(root)),'status':status},indent=2)+'\n')
    print(json.dumps({'status':status,'run':str(run)}))
    return 0 if status=='passed' else 1


if __name__=='__main__':
    raise SystemExit(main())
