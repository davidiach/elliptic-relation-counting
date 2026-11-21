"""Pinned dependencies and deterministic research records."""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parent
VENDOR = ROOT / 'vendor'
sys.path.insert(0, str(VENDOR))
if not __debug__:
    raise RuntimeError('These mathematical diagnostics require Python assertions.')


def verify_sources():
    manifest = json.loads((ROOT / 'source_manifest.json').read_text())
    for relative, entry in manifest['files'].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != entry['sha256']:
            raise RuntimeError('Pinned source changed: ' + relative)
    return {name: entry['sha256'] for name, entry in manifest['files'].items()}


def save(name, result):
    result['source_hashes'] = verify_sources()
    destination = Path(os.environ.get('HILLCLIMB_RESULTS', str(ROOT / 'results'))).resolve()
    if not destination.is_relative_to(ROOT):
        raise RuntimeError('Result directory must remain within this research iteration.')
    destination.mkdir(parents=True, exist_ok=True)
    (destination / name).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k in ('status', 'counts', 'scope')}, indent=2))
