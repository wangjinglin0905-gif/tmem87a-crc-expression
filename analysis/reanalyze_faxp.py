"""Recompute descriptive carcinoma-normal FAXP contrasts from PXD041412."""
import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd

TARGETS = {
    "TMEM87A": "Q8NBN3",
    "TMEM63A": "O94886",
    "EPCAM": "P16422",
    "KRT18": "P05783",
    "PTPRC": "P08575",
    "COL1A1": "P02452",
}


def read_design(path):
    # PRIDE deposits this design as an OOXML workbook with a .tsv suffix.
    with zipfile.ZipFile(path) as archive:
        strings = []
        for item in ET.fromstring(archive.read("xl/sharedStrings.xml")):
            strings.append("".join(node.text or "" for node in item.findall(".//{*}t")))
        rows = ET.fromstring(archive.read("xl/worksheets/sheet1.xml")).find(".//{*}sheetData")
        records = []
        for row in rows:
            values = {}
            for cell in row:
                value = cell.findtext("{*}v")
                if value is not None:
                    values[cell.attrib["r"].rstrip("0123456789")] = (
                        strings[int(value)] if cell.get("t") == "s" else value
                    )
            if "A" in values:
                records.append([values.get(col) for col in "ABCDEF"])
    design = pd.DataFrame(records[1:], columns=records[0])
    assert list(design.columns) == ["Batch", "ID", "Patient", "Slide", "Region", "Rep"]
    assert design.Batch.is_unique and design.ID.is_unique
    design["ms_batch"] = design.Batch.str.split("_").str[0]
    design["group"] = design.Region.replace({"CC": "C", "PC": "C"})
    return design


def contrast(regions, levels, group_col="group", carcinoma="C"):
    subset = regions[regions[group_col].isin([carcinoma, "N"])]
    both = subset.groupby([*levels, group_col]).agg(
        n=("log2_intensity", "size"), n_quantified=("log2_intensity", "count"),
        median_log2=("log2_intensity", "median"),
    ).unstack(group_col)
    both.columns = [f"{stat}_{group}" for stat, group in both.columns]
    both = both.reset_index()
    both["delta_log2"] = both[f"median_log2_{carcinoma}"] - both["median_log2_N"]
    return both


def main(source, output):
    design = read_design(source / "20230404_FAXP_CRC_Batch_design.tsv")
    matrix = pd.read_csv(source / "20230404_CRC_report.pg_matrix_DIA.tsv", sep="\t", low_memory=False)
    assay_cols = matrix.columns[5:]
    run_keys = {
        column: re.search(r"_(b\d+_(?:\d+|pool))_", column.rsplit("\\", 1)[-1])[1]
        for column in assay_cols
    }
    assert len(set(run_keys.values())) == len(assay_cols)
    map_report = {
        "matrix_protein_groups": len(matrix),
        "matrix_runs": len(run_keys),
        "design_runs": len(design),
        "design_missing_in_matrix": sorted(set(design.Batch) - set(run_keys.values())),
        "matrix_missing_in_design": sorted(set(run_keys.values()) - set(design.Batch)),
    }
    assert set(map_report["design_missing_in_matrix"]) == {"b1_1"}
    assert set(map_report["matrix_missing_in_design"]) == {f"b{n}_pool" for n in range(1, 7)}

    measures = []
    for gene, accession in TARGETS.items():
        row = matrix[
            matrix.Genes.astype(str).str.split(";").map(lambda genes: gene in genes)
        ]
        assert len(row) == 1 and row.iloc[0]["Protein.Group"] == accession
        assert row.iloc[0]["Protein.Ids"] == accession and row.iloc[0]["Genes"] == gene
        values = pd.to_numeric(row.iloc[0][assay_cols], errors="raise").astype(float).rename(index=run_keys)
        assert not np.isinf(values).any() and not (values < 0).any()
        sample = design.assign(
            gene=gene, accession=accession, intensity=design.Batch.map(values),
            in_matrix=design.Batch.isin(values.index),
        )
        sample["log2_intensity"] = np.log2(sample.intensity.where(sample.intensity > 0))
        measures.append(sample)
    measures = pd.concat(measures, ignore_index=True)
    observed = measures[measures.in_matrix].copy()
    levels = ["gene", "Patient", "Slide"]
    slide = contrast(observed, levels)
    batch = contrast(observed, [*levels, "ms_batch"])
    patient = (
        slide.dropna(subset=["delta_log2"]).groupby(["gene", "Patient"], as_index=False)
        .agg(slides=("Slide", "size"), delta_log2=("delta_log2", "mean"))
    )
    subtype = pd.concat(
        [
            contrast(observed[observed.Patient.eq("P1")], levels, "Region", key)
            .rename(columns={f"{stat}_{key}": f"{stat}_C" for stat in ["n", "n_quantified", "median_log2"]})
            .assign(carcinoma=key)
            for key in ["CC", "PC"]
        ],
        ignore_index=True,
    )
    dda = pd.read_csv(source / "20230404_CRC_combined_protein_DDA.tsv", sep="\t", low_memory=False)
    dda = dda.loc[dda["Gene"].isin(["TMEM87A", "TMEM63A"]), [
        "Protein ID", "Gene", "Combined Total Peptides",
        "Combined Unique Spectral Count", "Combined Total Spectral Count",
    ]]
    assert len(dda) == 2 and dda["Gene"].map(TARGETS).eq(dda["Protein ID"]).all()
    coverage = observed.groupby(["gene", "Patient", "Region"], as_index=False).agg(
        n_assayed=("intensity", "size"), n_quantified=("log2_intensity", "count"),
    )
    map_report["assayed_tissue_punches"] = int(observed.Batch.nunique())
    map_report["patients_with_carcinoma_and_normal"] = sorted(patient.Patient.unique())

    output.mkdir(exist_ok=True, parents=True)
    for name, frame in {
        "faxp_spot_values": measures, "faxp_coverage": coverage,
        "faxp_by_slide": slide, "faxp_by_patient": patient,
        "faxp_same_batch_sensitivity": batch, "faxp_carcinoma_subtypes": subtype,
        "faxp_dda_identification": dda,
    }.items():
        frame.to_csv(output / f"{name}.tsv", sep="\t", index=False, float_format="%.9g")
    (output / "faxp_join_qa.json").write_text(json.dumps(map_report, indent=2) + "\n", encoding="utf-8")
    print(patient.to_string(index=False))
    print(slide[slide.gene.isin(["TMEM87A", "TMEM63A"])].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.source, args.output)
