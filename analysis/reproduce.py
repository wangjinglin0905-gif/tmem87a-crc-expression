"""Reproduce core contrasts from the three extracted Additional files."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, spearmanr, wilcoxon


GENES = ("TMEM87A", "TMEM63A")
PAIRED_METRICS = ("log2_mean_ratio", "difference_cp10k", "detection_difference_pp")


def read(folder, name):
    return pd.read_csv(folder / name, sep="\t")


def paired_check(patients, summary, keys):
    checks = []
    for group, pairs in patients.groupby(keys):
        for metric in PAIRED_METRICS:
            values = pairs[metric].dropna()
            checks.append(dict(zip(keys, group), metric=metric, effect_check=values.mean(), n_check=len(values),
                               p_check=wilcoxon(values, method="auto").pvalue))
    checks = pd.DataFrame(checks)
    for _, subset in checks.groupby(keys[1:] + ["metric"]):
        checks.loc[subset.index, "q_check"] = false_discovery_control(subset.p_check.to_numpy())
    summary = summary[summary.metric.isin(PAIRED_METRICS)]
    matched = checks.merge(summary, on=keys + ["metric"], validate="one_to_one")
    assert len(matched) == len(checks) == len(summary)
    for name in ("n", "effect", "p", "q"):
        np.testing.assert_allclose(matched[name + "_check"], matched[name], atol=1e-12)
    return matched


def overlap(data, group, low, high, minimum):
    data = data[data[group].isin([low, high])].copy()
    cuts = np.unique(np.quantile(data.library, np.linspace(0, 1, 11)))
    data["depth"] = pd.cut(data.library, cuts, include_lowest=True, labels=False)
    counts = data.groupby(["depth", group]).size().unstack(fill_value=0).reindex(columns=[low, high], fill_value=0)
    counts = counts[counts.min(axis=1) >= minimum]
    eligible = data[data.depth.isin(counts.index)]
    assert len(counts), "No common support"
    for gene in GENES:
        profile = eligible.assign(detection=eligible[gene].gt(0)).groupby(["depth", group])[[gene, "detection"]].mean().unstack()
        weights = counts.min(axis=1)
        yield dict(gene=gene, n_strata=len(weights), n_common_weight=weights.sum(),
                   n_low=counts[low].sum(), n_high=counts[high].sum(),
                   n_low_original=data[group].eq(low).sum(), n_high_original=data[group].eq(high).sum(),
                   difference_cp10k=np.average(profile[gene][high] - profile[gene][low], weights=weights),
                   detection_difference_pp=100 * np.average(profile.detection[high] - profile.detection[low], weights=weights))


def paired_spatial_panel(spatial, depth, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8.1, 3.8))
    fig.subplots_adjust(left=.10, right=.98, top=.95, bottom=.29, wspace=.39)
    sections = sorted(spatial.section.unique())
    colors = dict(zip(GENES, ("#2468A2", "#BD5925")))
    for ax, (metric, label) in zip(axes, (("difference_cp10k", "Difference in mean CP10K"),
                                         ("detection_difference_pp", "Detection difference\n(percentage points)"))):
        ax.axhline(0, ls=":", lw=.8, color="#777777")
        for i, (gene, color) in enumerate(colors.items()):
            x = np.arange(1, len(sections) + 1) + (-.15 if i == 0 else .15)
            first = spatial.query("gene == @gene and view == 'all'").set_index("section").loc[sections, metric]
            units = [section.split("_")[1] for section in sections]
            second = depth.query("gene == @gene and contrast == 'spatial_proxy'").set_index("unit").loc[units, metric]
            ax.plot(np.vstack([x - .046, x + .046]), np.vstack([first, second]), color=color, lw=.9, alpha=.65)
            ax.scatter(x - .046, first, color=color, marker="o", s=20, label=gene + " · quartiles")
            ax.scatter(x + .046, second, edgecolor=color, facecolor="white", marker="o", s=20, label=gene + " · depth overlap")
        ax.set(xticks=range(1, 5), xlabel="Section (patient)", ylabel=label)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, ncol=2, loc="lower center", bbox_to_anchor=(.54, .012))
    fig.savefig(output / "GSE226997_spatial_contrasts.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="Parent directory of the three extracted Additional_file_* directories")
    parser.add_argument("--output", type=Path, required=True, help="A writable directory for regenerated results")
    args = parser.parse_args()
    f1, f2, f3 = (args.data / f"Additional_file_{n}" for n in (1, 2, 3))
    args.output.mkdir(parents=True, exist_ok=True)
    paired = paired_check(read(f2, "GSE166555_patient_pairs.tsv").query("cutoff==20 and weighting=='pooled_cells'"),
                          read(f2, "GSE166555_paired_summary.tsv").query("cutoff==20 and weighting=='pooled_cells'"), ["gene"])
    cna = paired_check(read(f2, "GSE166555_SCNA_patient_pairs.tsv").query("cutoff==20"),
                       read(f2, "GSE166555_SCNA_paired_summary.tsv").query("cutoff==20"), ["gene", "contrast"])
    pd.concat([paired.assign(contrast="Tumour vs adjacent"), cna], ignore_index=True).to_csv(args.output / "paired_recomputed.tsv", sep="\t", index=False)

    cptac_pairs = read(f1, "cptac_patient_targets.tsv")
    cptac_expected = read(f1, "cptac_paired_protein.tsv").set_index("gene")
    cptac_rna = read(f1, "cptac_rna_protein.tsv").set_index("gene")
    for gene, values in cptac_pairs.groupby("gene"):
        pairs = values[["tumour", "normal"]].dropna()
        rna = values[["rna", "tumour"]].dropna()
        assert len(pairs) == cptac_expected.at[gene, "n_complete_pairs"]
        np.testing.assert_allclose((pairs.tumour - pairs.normal).mean(), cptac_expected.at[gene, "mean_difference"])
        np.testing.assert_allclose(spearmanr(rna.rna, rna.tumour).statistic, cptac_rna.at[gene, "rho"])

    depth_rows = []
    epi = read(f2, "GSE166555_epithelial_targets.tsv.gz")
    for label, frame in [("all_epithelial", epi), ("CNA_vs_adjacent", epi[epi.sample_origin.eq("Normal") | epi.cna_clone.eq("CNA")])]:
        for patient, rows in frame.groupby("case_id"):
            if rows.sample_origin.value_counts().reindex(["Normal", "Tumor"], fill_value=0).min() >= 20:
                depth_rows.extend(dict(unit=patient, contrast=label, **row) for row in overlap(rows, "sample_origin", "Normal", "Tumor", 5))
    spatial = read(f3, "GSE226997_spatial_proxy_summary.tsv")
    for path in sorted(f3.glob("GSE226997_*_spot_targets.tsv.gz")):
        full = read(f3, path.name)
        section = path.name.removesuffix("_spot_targets.tsv.gz")
        unit = section.split("_")[1]
        for view, data in (("all", full), ("middle_80pct_library", full[full.depth_middle_80])):
            means = data.groupby("stratum")[list(GENES)].mean()
            detection = data[list(GENES)].gt(0).groupby(data.stratum).mean()
            expect = spatial.query("section == @section and view == @view").set_index("gene").loc[list(GENES)]
            np.testing.assert_allclose(means.loc["high"] - means.loc["low"], expect.difference_cp10k, atol=1e-12)
            np.testing.assert_allclose(100 * (detection.loc["high"] - detection.loc["low"]), expect.detection_difference_pp, atol=1e-12)
        depth_rows.extend(dict(unit=unit, contrast="spatial_proxy", **row) for row in overlap(full, "stratum", "low", "high", 10))
    depth = pd.DataFrame(depth_rows)
    expected = read(f3, "depth_overlap_diagnostics.tsv")
    index = ["unit", "contrast", "gene"]
    pd.testing.assert_frame_equal(depth.set_index(index).sort_index(), expected.set_index(index).sort_index()[depth.columns[3:]], check_exact=False, atol=1e-11)
    depth.to_csv(args.output / "depth_recomputed.tsv", sep="\t", index=False)
    paired_spatial_panel(spatial, depth, args.output)
    receipt = dict(paired_estimates_checked=len(paired)+len(cna), cptac_proteins=2, spatial_donors=4,
                   depth_gene_patient_comparisons=len(depth), scope="Derived target and patient data")
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
