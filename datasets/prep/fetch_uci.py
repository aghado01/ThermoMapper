#!/usr/bin/env python3
"""
fetch_uci.py — provenance + regenerator for the UCI benchmark CSVs in this
directory (landsat.csv, isolet.csv.gz). Run from the repo root:

    python datasets/prep/fetch_uci.py

Network note (observed in the portable env): `archive.ics.uci.edu` resolves but
`openml.org` does NOT (selective DNS). UCI ships ISOLET as Unix-`compress` (.Z)
members inside the zip, which Python's stdlib zipfile cannot inflate — so we
shell out to `gzip -d` (handles .Z). The datasets themselves are committed; this
script exists for provenance and reproducibility, not as a build step.

Sources:
  Landsat (Statlog Landsat Satellite) — UCI id 146
  ISOLET (Isolated Letter Speech Recognition) — UCI id 54
"""
import io
import os
import ssl
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

LANDSAT_URL = "https://archive.ics.uci.edu/static/public/146/statlog+landsat+satellite.zip"
ISOLET_URL = "https://archive.ics.uci.edu/static/public/54/isolet.zip"

HERE = os.path.dirname(os.path.abspath(__file__))
DATASETS = os.path.dirname(HERE)  # datasets/
_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE


def _download(url: str) -> bytes:
    return urllib.request.urlopen(url, timeout=180, context=_CTX).read()


def build_landsat() -> None:
    """sat.trn + sat.tst → landsat.csv (space→comma; 36 int features + int label
    {1,2,3,4,5,7}; class 6 'mixture' is absent in the source). Full 6435 rows."""
    z = zipfile.ZipFile(io.BytesIO(_download(LANDSAT_URL)))
    rows = []
    for member in ("sat.trn", "sat.tst"):
        for line in z.read(member).decode().splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            assert len(parts) == 37, len(parts)
            rows.append(",".join(parts))
    out = os.path.join(DATASETS, "landsat.csv")
    with open(out, "w", newline="") as f:
        f.write("\n".join(rows) + "\n")
    print(f"landsat.csv: {len(rows)} rows")


def build_isolet() -> None:
    """isolet1+2+3+4.data + isolet5.data → isolet.csv.gz (full 7797 rows; 617
    continuous features in [-1,1] + int letter label 1..26, parsed from the
    source's float form e.g. '1.'). Gzipped (~10MB vs 35MB raw)."""
    z = zipfile.ZipFile(io.BytesIO(_download(ISOLET_URL)))
    with tempfile.TemporaryDirectory() as work:
        data_files = []
        for member in ("isolet1+2+3+4.data.Z", "isolet5.data.Z"):
            zpath = os.path.join(work, os.path.basename(member))
            with open(zpath, "wb") as f:
                f.write(z.read(member))
            # .Z is Unix compress (LZW); gzip -d handles it, stdlib cannot.
            subprocess.run(["gzip", "-d", "-f", zpath], check=True)
            data_files.append(zpath[:-2])  # strip .Z

        rows = []
        for df in data_files:
            for line in open(df):
                line = line.strip()
                if not line:
                    continue
                parts = [p.strip() for p in line.split(",")]
                assert len(parts) == 618, len(parts)
                label = int(round(float(parts[-1])))  # '1.' -> 1
                rows.append(",".join(parts[:-1]) + "," + str(label))

    csv_path = os.path.join(DATASETS, "isolet.csv")
    with open(csv_path, "w", newline="") as f:
        f.write("\n".join(rows) + "\n")
    subprocess.run(["gzip", "-9", "-f", csv_path], check=True)  # -> isolet.csv.gz
    print(f"isolet.csv.gz: {len(rows)} rows")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "landsat"):
        build_landsat()
    if which in ("all", "isolet"):
        build_isolet()
