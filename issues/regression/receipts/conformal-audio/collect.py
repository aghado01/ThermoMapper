"""Retain compact receipts with actual-output and LF-copy hashes."""
import argparse
import hashlib
import json
import pathlib


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


p = argparse.ArgumentParser()
p.add_argument("--run", type=pathlib.Path, required=True)
p.add_argument("--output", type=pathlib.Path, required=True)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
records = []
local_only = {"audio-tests.jsonl", "audio-predictions.jsonl", "risk-conversations.jsonl"}
for source in sorted((a.run / "results").glob("*")) + sorted((a.run / "analysis").glob("*.json")):
    if not source.is_file():
        continue
    record = {"actual": source.as_posix(), "bytes": source.stat().st_size, "sha256": sha(source)}
    if source.name not in local_only and not source.name.startswith("training-"):
        target = a.output / source.name
        target.write_bytes(source.read_bytes().replace(b"\r\n", b"\n"))
        record.update(receipt=target.as_posix(), receipt_sha256_lf=sha(target))
    else:
        record["raw_local_only"] = True
    records.append(record)
(a.output / "verification.json").write_bytes((json.dumps({"files": records}, indent=2) + "\n").encode())
print(json.dumps({"files": len(records), "tracked_copies": sum("receipt" in r for r in records)}))
