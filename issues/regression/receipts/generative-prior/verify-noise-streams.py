"""Compare G5's Python integer streams with an isolated compiled C# assembly."""
import argparse
import json
import math
import pathlib
import runpy

p = argparse.ArgumentParser()
p.add_argument("--reference", type=pathlib.Path, required=True)
p.add_argument("--output", type=pathlib.Path, required=True)
a = p.parse_args()
reference = json.loads(a.reference.read_text(encoding="utf-8-sig"))
code = runpy.run_path(str(pathlib.Path(__file__).with_name("g5-noise.py")))
for case in reference["cases"]:
    assert code["derive"](case["seed"], 10) == case["children"]
    rng = code["Xoshiro"](case["seed"])
    assert [str(rng.next()) for _ in range(10)] == case["words"]
    rng = code["Xoshiro"](case["seed"])
    assert [rng.below(30) for _ in range(10)] == case["bounded"]
rows = [{"k": k, "m": k * ratio, "noise": noise,
         "rmse": math.exp(0.5 * (k / 3 + 1.25 * math.log(1 / ratio) + 0.75 * math.log(noise)))}
        for k in (1, 2, 3) for ratio in (2, 4, 8, 16) for noise in (0.01, 0.04, 0.16)]
fit = code["slopes"](rows)
assert abs(fit["ratio_slope"] - 1.25) < 1e-12 and abs(fit["noise_slope"] - 0.75) < 1e-12
single_noise = code["slopes"]([r for r in rows if r["noise"] == 0.04])
assert abs(single_noise["ratio_slope"] - 1.25) < 1e-12 and single_noise["noise_slope"] is None
result = {"passed": True, "cases": len(reference["cases"]), "reference": str(a.reference),
          "known_regression_slopes": fit, "single_noise_control": single_noise,
          "assembly_sha256": reference["assembly_sha256"]}
a.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))
