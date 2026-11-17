"""Check portable G5 replay and reconstruct its stored per-cell readings."""
import argparse
import collections
import csv
import hashlib
import json
import pathlib
import statistics

p = argparse.ArgumentParser()
p.add_argument("--primary", type=pathlib.Path, required=True)
p.add_argument("--replay", type=pathlib.Path, required=True)
p.add_argument("--output", type=pathlib.Path, required=True)
a = p.parse_args()
primary = json.loads((a.primary / "summary.json").read_text())
replay = json.loads((a.replay / "summary.json").read_text())
assert primary["stages"] == replay["stages"]
assert (a.primary / "oracle-rows.csv").read_bytes() == (a.replay / "oracle-rows.csv").read_bytes()
with (a.primary / "oracle-rows.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 1960
clusters = collections.defaultdict(set)
for r in rows:
    clusters[r["stage"], int(r["k"])].add(int(r["rep"]))
for stage in primary["stages"]:
    for cell in stage["cells"]:
        selected = [r for r in rows if (r["stage"], int(r["k"]), float(r["noise"]), int(r["ratio"])) ==
                    (stage["stage"], cell["k"], cell["noise"], cell["ratio"])]
        errors = [float(r["rmse"]) ** 2 for r in selected]
        assert len(errors) == cell["rows"]
        assert statistics.median(errors) == cell["median_mse"]
        assert cell["converged_rows"] == sum(int(r["converged"]) == int(r["runs"]) for r in selected)
        assert cell["rounding_lower_median"] <= cell["median_mse"] <= cell["rounding_upper_median"]
    converged = [r for r in rows if r["stage"] == stage["stage"] and int(r["converged"]) == int(r["runs"])]
    assert len(converged) == stage["converged_rows"]
result = {"passed": True, "exact_stage_payload_replay": True, "exact_portable_rows_replay": True,
          "portable_sha256": hashlib.sha256((a.primary / "oracle-rows.csv").read_bytes()).hexdigest(),
          "clusters_per_stage_k": [{"stage": s, "k": k, "clusters": len(v)} for (s, k), v in clusters.items()],
          "primary_settings_sha256": primary["settings_sha256"], "replay_settings_sha256": replay["settings_sha256"],
          "primary_source_sha256": primary["analysis_source_sha256"], "replay_source_sha256": replay["analysis_source_sha256"]}
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))
