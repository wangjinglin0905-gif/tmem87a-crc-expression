"""Reproduce the post hoc patient-level RNA decomposition from Additional files 1 and 2."""
import argparse
import gzip
import json
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


def table(folder, number, name):
    subdir = f"Additional_file_{number}"
    with ZipFile(folder / (subdir + ".zip")) as archive:
        content = archive.read(f"{subdir}/{name}")
    if name.endswith(".gz"):
        content = gzip.decompress(content)
    return pd.read_csv(BytesIO(content), sep="\t")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("supplements", type=Path, help="Folder containing the two numbered ZIP files")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    discovery = table(args.supplements, 1, "epithelial_patient_pairs.tsv")
    discovery = discovery.query("view == 'all_processing' and min_epithelial_cells == 20")
    external = table(args.supplements, 2, "GSE166555_patient_pairs.tsv")
    external = external.query("weighting == 'pooled_cells' and cutoff == 20").assign(dataset="GSE166555")
    keys = ["dataset", "gene", "patient"]
    cols = [s + "_" + m for s in ("tumour", "normal") for m in ("n", "mean", "detection")]
    pairs = pd.concat([discovery[keys + cols], external[keys + cols]]).sort_values(keys).reset_index(drop=True)
    assert not pairs.duplicated(keys).any()
    assert set(pairs.gene) == {"TMEM63A", "TMEM87A"}
    assert (pairs[["tumour_n", "normal_n"]] >= 20).all().all()
    for side in ("tumour", "normal"):
        p = pairs[side + "_detection"]
        x = pairs[side + "_mean"]
        assert p.between(0, 1).all() and x.ge(0).all() and p.eq(0).equals(x.eq(0))
        pairs[side + "_mean_detected"] = x / p.where(p.gt(0))
    for name, suffix in [("total", "mean"), ("detection", "detection"), ("detected_cell_mean", "mean_detected")]:
        a, b = (pairs[s + "_" + suffix] for s in ("tumour", "normal"))
        pairs[name] = np.log2(a.where(a.gt(0)) / b.where(b.gt(0)))
    metrics = ["total", "detection", "detected_cell_mean"]
    valid = pairs.dropna(subset=metrics)
    residual = abs(valid.total - valid.detection - valid.detected_cell_mean).max()
    assert residual < 1e-10
    cells = table(args.supplements, 2, "GSE166555_epithelial_targets.tsv.gz")
    direct = 0
    for record in external.itertuples():
        patient = cells[cells.case_id == record.patient]
        for side, origin in [("tumour", "Tumor"), ("normal", "Normal")]:
            values = patient.loc[patient.sample_origin.eq(origin), record.gene]
            assert len(values) == getattr(record, side + "_n")
            x, p = getattr(record, side + "_mean"), getattr(record, side + "_detection")
            np.testing.assert_allclose([values.mean(), values.gt(0).mean(), values[values.gt(0)].mean()],
                                       [x, p, x / p if p else np.nan], rtol=1e-9, atol=1e-12, equal_nan=True)
            direct += 1
    results = []
    for (dataset, gene), group in pairs.groupby(keys[:2], sort=True):
        data = group[metrics].dropna().to_numpy()
        # Resample the same patients for all three components. Intervals are marginal.
        rng = np.random.default_rng(20261125)
        sampled = data[rng.integers(0, len(data), size=(5000, len(data)))].mean(axis=1)
        for col, metric in enumerate(metrics):
            low, high = np.quantile(sampled[:, col], [0.025, 0.975])
            results.append(dict(dataset=dataset, gene=gene, metric=metric, eligible=len(group), defined=len(data),
                                mean=data[:, col].mean(), ci_low=low, ci_high=high))
    args.output.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.output / "patient_components.tsv", sep="\t", index=False, float_format="%.12g")
    pd.DataFrame(results).to_csv(args.output / "cohort_components.tsv", sep="\t", index=False, float_format="%.12g")
    print(json.dumps(dict(patient_gene_pairs=len(pairs), defined=len(valid), max_identity_error=float(residual),
                          external_patient_gene_tissue_checks=direct)))


if __name__ == "__main__":
    main()
