"""Check exact historical archive members without executing them."""
from pathlib import Path
import ast
import hashlib
import json


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    archive_paths = {
        path.relative_to(root).as_posix()
        for path in root.glob("[0-9][0-9]_*/**/*")
        if path.is_file()
    }
    manifest_paths = set()
    scripts = 0
    for cp in manifest["checkpoints"]:
        for record in cp["files"]:
            member = record["path"]
            if member in manifest_paths:
                raise SystemExit(f"Duplicate manifest path: {member}")
            if member not in archive_paths:
                raise SystemExit(f"Missing archive member: {member}")
            manifest_paths.add(member)
            p = root / member
            data = p.read_bytes()
            if hashlib.sha256(data).hexdigest() != record["sha256"] or len(data) != record["bytes"]:
                raise SystemExit(f"Archive integrity failure: {p}")
            blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            if blob != record["git_blob"]:
                raise SystemExit(f"Git blob mismatch: {p}")
            if p.suffix == ".py":
                ast.parse(data.decode("utf-8"), filename=str(p))
                scripts += 1
    unlisted = archive_paths - manifest_paths
    if unlisted:
        raise SystemExit("Archive members missing from manifest: " + ", ".join(sorted(unlisted)))
    print(f"Verified {len(manifest['checkpoints'])} checkpoints, {len(manifest_paths)} archived files, {scripts} Python source files; no historical scripts executed.")


if __name__ == "__main__":
    main()
