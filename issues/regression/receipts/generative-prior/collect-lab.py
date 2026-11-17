"""Copy compact lab receipts, recording actual outputs and exact replay checks."""
import argparse
import hashlib
import json
import pathlib
import shutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True, type=pathlib.Path)
    p.add_argument("--output", required=True, type=pathlib.Path)
    p.add_argument("--replay", type=pathlib.Path)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    records = []
    for source in sorted((a.run / "results").glob("*")):
        if not source.is_file():
            continue
        record = {"actual": source.as_posix(), "bytes": source.stat().st_size, "sha256": sha(source)}
        if source.name not in ("posterior-draws.jsonl", "alignment-states.jsonl"):
            target = a.output / source.name
            target.write_bytes(source.read_bytes().replace(b"\r\n", b"\n"))
            record["receipt"] = target.as_posix()
            record["receipt_sha256_lf"] = sha(target)
        else:
            record["raw_states_local_only"] = True
        if a.replay:
            replay = a.replay / "results" / source.name
            record["replay"] = replay.as_posix()
            record["replay_sha256"] = sha(replay)
            record["exact_replay"] = record["sha256"] == record["replay_sha256"]
            assert record["exact_replay"], source
        records.append(record)
    for source in sorted((a.run / "analysis").glob("*.json")):
        (a.output / source.name).write_bytes(source.read_bytes().replace(b"\r\n", b"\n"))
        records.append({"actual": source.as_posix(), "receipt": (a.output / source.name).as_posix(),
                        "bytes": source.stat().st_size, "sha256": sha(source),
                        "receipt_sha256_lf": sha(a.output / source.name)})
    (a.output / "verification.json").write_bytes((json.dumps({"files": records}, indent=2) + "\n").encode())
    print(json.dumps({"files": len(records), "copied_bytes": sum(r["bytes"] for r in records if "receipt" in r),
                      "exact_replay": bool(a.replay)}))


if __name__ == "__main__":
    main()
