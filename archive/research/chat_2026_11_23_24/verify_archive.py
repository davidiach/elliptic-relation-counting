"""Check exact historical archive members without executing them."""
from pathlib import Path
import ast
import hashlib
import json


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    total = scripts = 0
    for cp in manifest["checkpoints"]:
        for record in cp["files"]:
            p = root / record["path"]
            data = p.read_bytes()
            if hashlib.sha256(data).hexdigest() != record["sha256"] or len(data) != record["bytes"]:
                raise SystemExit(f"Archive integrity failure: {p}")
            blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            if blob != record["git_blob"]:
                raise SystemExit(f"Git blob mismatch: {p}")
            if p.suffix == ".py":
                ast.parse(data.decode("utf-8"), filename=str(p))
                scripts += 1
            total += 1
    print(f"Verified {len(manifest['checkpoints'])} checkpoints, {total} archived files, {scripts} Python source files; no historical scripts executed.")


if __name__ == "__main__":
    main()
