"""Resolve example paths without downloading references or running the design."""

import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    params = json.loads((root / "params.template.json").read_text())
    for key in ("fg_genomes", "bg_genomes", "fg_prefixes", "bg_prefixes"):
        params[key] = [str(root / value) for value in params[key]]
    params["data_dir"] = str(root / params["data_dir"])
    for value in params["fg_genomes"] + params["bg_genomes"]:
        if not Path(value).is_file():
            raise FileNotFoundError(f"Missing reference: {value}. See README.md.")
    (root / "work").mkdir(exist_ok=True)
    destination = root / "params.json"
    destination.write_text(json.dumps(params, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
