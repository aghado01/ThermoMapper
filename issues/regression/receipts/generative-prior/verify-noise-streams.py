"""Compare G5's Python integer streams with an isolated compiled C# assembly."""
import argparse
import json
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
result = {"passed": True, "cases": len(reference["cases"]), "reference": str(a.reference),
          "assembly_sha256": reference["assembly_sha256"]}
a.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))
