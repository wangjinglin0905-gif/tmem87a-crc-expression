"""Verify archive integrity and reproduce selected results from the numbered supplements."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path, PurePosixPath
from zipfile import ZipFile


def verify(data, expected, name):
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"SHA-256 mismatch for {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="An empty output directory")
    output = parser.parse_args().output.resolve()
    release = Path(__file__).resolve().parent
    for name, value in json.loads((release / "sha256.json").read_text(encoding="utf-8")).items():
        verify((release / name).read_bytes(), value, name)
    if output.exists() and any(output.iterdir()):
        parser.error("Select an empty directory to preserve existing results.")
    output.mkdir(parents=True, exist_ok=True)
    for n in range(1, 6):
        label = f"Additional_file_{n}"
        with ZipFile(release / "supplement" / (label + ".zip")) as archive:
            manifest = json.loads(archive.read(f"{label}/sha256.json"))
            expected_names = {f"{label}/{k}" for k in manifest} | {f"{label}/sha256.json"}
            if set(archive.namelist()) != expected_names or len(archive.namelist()) != len(expected_names):
                raise ValueError(f"Unexpected archive members in {label}")
            for name in archive.namelist():
                path = PurePosixPath(name)
                if len(path.parts) != 2 or path.parts[0] != label or path.name in (".","..") or "\\" in name or ":" in name:
                    raise ValueError(f"Invalid archive member: {name}")
                data = archive.read(name)
                if path.name != "sha256.json":
                    verify(data, manifest[path.name], name)
                dest = output / label / path.name
                dest.parent.mkdir(exist_ok=True)
                dest.write_bytes(data)
    commands = [
        [output / "Additional_file_3/reproduce.py", "--data", output, "--output", output / "core"],
        [output / "Additional_file_4/verify_proteomics.py"],
        [output / "Additional_file_5/decompose_expression.py", release / "supplement", output / "decomposition"],
    ]
    for command in commands:
        subprocess.run([sys.executable, *map(str, command)], cwd=output, check=True)
    for name in ("cohort_components.tsv", "patient_components.tsv"):
        actual = (output / "decomposition" / name).read_text(encoding="utf-8")
        expected = (output / "Additional_file_5" / name).read_text(encoding="utf-8")
        if actual != expected:
            raise ValueError(f"Regenerated table differs: {name}")
    print("Verified all five supplement archives and reproduced the selected results.")


if __name__ == "__main__":
    main()
