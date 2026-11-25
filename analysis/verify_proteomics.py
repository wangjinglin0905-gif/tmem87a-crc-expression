"""Recalculate the external proteomics summaries from the bundled measurements."""
from hashlib import file_digest
from itertools import product
import json
from math import fsum, isclose, log2, sqrt
from pathlib import Path
from statistics import median, stdev

import pandas as pd
from scipy.stats import false_discovery_control, t


def close(actual, expected):
    assert isclose(actual, expected, rel_tol=2e-7, abs_tol=1e-7), (actual, expected)


def main():
    root = Path(__file__).resolve().parent
    for name, expected in json.loads((root / "sha256.json").read_text(encoding="utf-8")).items():
        with (root / name).open("rb") as handle:
            assert file_digest(handle, "sha256").hexdigest() == expected, name

    lcm = pd.read_csv(root / "lcm_pairs.tsv", sep="\t")
    summaries = pd.read_csv(root / "lcm_estimates.tsv", sep="\t")
    for row in summaries.itertuples():
        pairs = lcm[lcm.gene.eq(row.gene)]
        if row.analysis == "exclude_case_10":
            pairs = pairs[pairs.case.ne(10)]
        delta = [log2(a / b) for a, b in zip(pairs.carcinoma_lfq, pairs.normal_lfq)]
        n = len(delta)
        assert n == row.n_complete
        assert sum(x < 0 for x in delta) == row.n_lower_in_tumor
        assert sum(x > 0 for x in delta) == row.n_higher_in_tumor
        mean = fsum(delta) / n
        close(mean, row.mean_log2_difference)
        close(median(delta), row.median_log2_difference)
        close(2**mean, row.geometric_mean_ratio)
        width = t.ppf(0.975, n - 1) * stdev(delta) / sqrt(n)
        close(mean - width, row.ci95_lower_t)
        close(mean + width, row.ci95_upper_t)
        flips = (fsum(x * sign for x, sign in zip(delta, signs)) for signs in product([-1, 1], repeat=n))
        close(sum(abs(value) >= abs(fsum(delta)) - 1e-10 for value in flips) / 2**n, row.exact_p_sign_flip)
    primary = summaries[summaries.bh_q_two_targets.notna()]
    assert set(primary.gene) == {"TMEM87A", "TMEM63A"} and primary.analysis.eq("all_pairs").all()
    for recomputed, stored in zip(false_discovery_control(primary.exact_p_sign_flip), primary.bh_q_two_targets):
        close(recomputed, stored)

    spot = pd.read_csv(root / "faxp_spot_values.tsv", sep="\t")
    measured = spot[spot.in_matrix & spot.intensity.gt(0)]
    patient = pd.read_csv(root / "faxp_by_patient.tsv", sep="\t")
    for row in patient.itertuples():
        data = measured[measured.gene.eq(row.gene) & measured.Patient.eq(row.Patient)]
        values = []
        for _, section in data.groupby("Slide"):
            a = section.loc[section.group.eq("C"), "intensity"]
            b = section.loc[section.group.eq("N"), "intensity"]
            if len(a) and len(b):
                values.append(median(map(log2, a)) - median(map(log2, b)))
        assert len(values) == row.slides
        close(fsum(values) / len(values), row.delta_log2)
    print(json.dumps({
        "archive_member_hashes": "verified",
        "LCM_summary_checks": len(summaries),
        "FAXP_patient_checks": len(patient),
        "LCM_target_BH_adjustments": len(primary),
    }))


if __name__ == "__main__":
    main()
