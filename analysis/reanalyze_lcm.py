"""Read MaxQuant LFQ from PXD002137 and summarize patient-matched contrasts."""
import argparse
import itertools
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd
from scipy import stats

from reanalyze_faxp import TARGETS


def summarize(delta):
    values = np.asarray(delta, dtype=float)
    n = len(values)
    mean = float(values.mean())
    half_width = float(stats.t.ppf(0.975, n - 1) * stats.sem(values))
    permutations = np.array(list(itertools.product([-1, 1], repeat=n)))
    null = np.abs(permutations @ values / n)
    return {
        "n_complete": n, "n_higher_in_tumor": int((values > 0).sum()),
        "n_lower_in_tumor": int((values < 0).sum()),
        "mean_log2_difference": mean, "geometric_mean_ratio": 2**mean,
        "median_log2_difference": float(np.median(values)),
        "ci95_lower_t": mean - half_width, "ci95_upper_t": mean + half_width,
        "exact_p_sign_flip": float(np.mean(null >= abs(mean) - 1e-12)),
    }


def main(source, output):
    with zipfile.ZipFile(source) as archive:
        assert archive.testzip() is None, "Corrupt ZIP member"
        with archive.open("experimentalDesignTemplate.txt") as handle:
            design = pd.read_csv(handle, sep="\t")
        with archive.open("proteinGroups.txt") as handle:
            groups = pd.read_csv(handle, sep="\t", low_memory=False, encoding="cp1252")

    # The source paper reports paired adjacent mucosa and interleaved C1/N1, C2/N2.
    # The deposit gives the complete CAN/NOR identifiers with the shared case number.
    matched = design[design.Experiment.str.fullmatch(r"[CN]\d+")].copy()
    details = matched.Name.str.extract(r"^(\d{8}).*_(CAN|NOR)(\d+)_")
    assert details.notna().all().all()
    assert (details[1].map({"CAN": "C", "NOR": "N"}).str.cat(details[2]) == matched.Experiment).all()
    matched["case"] = details[2].astype(int)
    matched["run_date"] = details[0]
    pairs = sorted(int(case) for case in matched["case"].unique())
    assert pairs == [1, 2, 3, 4, 5, 8, 9, 10]
    expected = {f"{group}{case}" for case in pairs for group in ["C", "N"]}
    assert set(matched.Experiment) == expected
    assert set(name.removeprefix("LFQ intensity ") for name in groups if re.match(
        r"^LFQ intensity [CN]\d+$", name,
    )) == expected

    flags = ["Reverse", "Contaminant", "Only identified by site"]
    clean = groups[groups[flags].fillna("").ne("+").all(axis=1)]
    row_info, pair_values, estimates = [], [], []
    for gene, accession in TARGETS.items():
        exact_names = clean["Gene names"].fillna("").str.split(";")
        candidates = clean[exact_names.map(lambda names: gene in names)]
        row = candidates[candidates["Protein IDs"].fillna("").str.contains(
            fr"(?:^|;){re.escape(accession)}(?:-|;|$)", regex=True
        )]
        assert len(row) == 1
        row = row.iloc[0]
        assert row["Gene names"] == gene and row["Unique peptides"] > 0
        row_info.append({
            "gene": gene, "accession": accession, "protein_ids": row["Protein IDs"],
            "majority_ids": row["Majority protein IDs"],
            "peptides": int(row["Peptides"]), "unique_peptides": int(row["Unique peptides"]),
        })
        records = []
        for case in pairs:
            tumor = float(row[f"LFQ intensity C{case}"])
            normal = float(row[f"LFQ intensity N{case}"])
            records.append({
                "gene": gene, "case": case, "carcinoma_lfq": tumor, "normal_lfq": normal,
                "carcinoma_peptides": int(row[f"Peptides C{case}"]),
                "normal_peptides": int(row[f"Peptides N{case}"]),
                "delta_log2": float(np.log2(tumor / normal)) if tumor > 0 and normal > 0 else np.nan,
            })
        pair_values.extend(records)
        for name, subset in [
            ("all_pairs", records), ("exclude_case_10", [item for item in records if item["case"] != 10])
        ]:
            delta = [item["delta_log2"] for item in subset if np.isfinite(item["delta_log2"])]
            estimates.append(dict(gene=gene, analysis=name, **summarize(delta)))

    estimates = pd.DataFrame(estimates)
    primary = estimates.analysis.eq("all_pairs") & estimates.gene.isin(["TMEM87A", "TMEM63A"])
    estimates.loc[primary, "bh_q_two_targets"] = stats.false_discovery_control(
        estimates.loc[primary, "exact_p_sign_flip"]
    )
    output.mkdir(exist_ok=True, parents=True)
    for title, frame in {
        "lcm_protein_groups": pd.DataFrame(row_info), "lcm_design": matched,
        "lcm_pairs": pd.DataFrame(pair_values), "lcm_estimates": estimates,
    }.items():
        frame.to_csv(output / f"{title}.tsv", sep="\t", index=False, float_format="%.9g")
    (output / "lcm_qa.json").write_text(json.dumps({
        "archive_crc_ok": True, "n_maxquant_groups": len(groups),
        "n_groups_after_flags": len(clean),
        "flag_columns": flags, "patient_pair_codes": pairs,
        "n_adenomas_excluded": int(design[design.Experiment.str.startswith("A")].Experiment.nunique()),
    }, indent=2) + "\n", encoding="utf-8")
    print(estimates[estimates.gene.isin(["TMEM87A", "TMEM63A"])].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.source, args.output)
