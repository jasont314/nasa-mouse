"""Split a skeletal-muscle expiMap AnnData query by OSDR material muscle group."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .analyze_muscle_targeted_modules import muscle_group
from .io import require_import


DEFAULT_INPUT = (
    "outputs/expimap_muscle_targeted_combined/input/"
    "osdr_skeletal_muscle_flt_gc_reactome_raw_counts.h5ad"
)
DEFAULT_OUTPUT_DIR = "outputs/expimap_muscle_targeted_combined/group_inputs"


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")

    adata = ad.read_h5ad(args.input)
    material_col = "study.characteristics.material type"
    if material_col not in adata.obs:
        raise SystemExit(f"Input missing `{material_col}` in obs.")
    adata.obs["muscle_group"] = adata.obs[material_col].map(muscle_group)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    outputs = []
    for group_name, obs in adata.obs.groupby("muscle_group", sort=True):
        sub = adata[obs.index.astype(str).tolist(), :].copy()
        counts = (
            sub.obs.groupby(["id.accession", "condition_inferred"], dropna=False)
            .size()
            .unstack(fill_value=0)
        )
        for condition in ["flight", "ground_control"]:
            if condition not in counts:
                counts[condition] = 0
        accessions_with_both = int(((counts["flight"] > 0) & (counts["ground_control"] > 0)).sum())
        n_flight = int((sub.obs["condition_inferred"] == "flight").sum())
        n_ground = int((sub.obs["condition_inferred"] == "ground_control").sum())
        out_path = output_dir / f"osdr_skeletal_muscle_{group_name}_flt_gc_reactome_plus_muscle_raw_counts.h5ad"
        if (
            n_flight >= args.min_per_condition
            and n_ground >= args.min_per_condition
            and accessions_with_both >= args.min_accessions
        ):
            sub.write_h5ad(out_path)
            outputs.append(str(out_path))
        rows.append(
            {
                "muscle_group": group_name,
                "samples": int(sub.n_obs),
                "flight": n_flight,
                "ground_control": n_ground,
                "accessions_with_both": accessions_with_both,
                "genes": int(sub.n_vars),
                "terms": int(len(sub.uns.get("terms", []))),
                "written": str(out_path) if str(out_path) in outputs else "",
            }
        )

    summary = pd.DataFrame(rows).sort_values(["accessions_with_both", "samples"], ascending=False)
    summary_path = output_dir / "muscle_group_input_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    manifest = {
        "input": str(args.input),
        "output_dir": str(output_dir),
        "min_per_condition": args.min_per_condition,
        "min_accessions": args.min_accessions,
        "groups_written": len(outputs),
        "outputs": outputs,
        "summary": str(summary_path),
    }
    manifest_path = output_dir / "muscle_group_input_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split skeletal-muscle expiMap input by muscle group.")
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-per-condition", type=int, default=3)
    parser.add_argument("--min-accessions", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
