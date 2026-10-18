import math
import pandas as pd
from pathlib import Path

INPUT = Path("data/gnomad_variant_frequencies.csv")
WEIGHTS = Path("data/census_scenario_weights.csv")
OUTPUT = Path("data/source_sensitivity_estimates.csv")

GENES = ["ABCG2", "CYP3A5", "CYP4F2", "NUDT15", "VKORC1"]
PANELS = {
    "all_five_panel": GENES,
    "reduced_panel": ["ABCG2", "CYP4F2", "NUDT15"],
}

def wilson(ac, an, z=1.96):
    p = ac / an
    z2 = z * z
    denominator = 1 + z2 / an
    centre = (p + z2 / (2 * an)) / denominator
    half = z * math.sqrt(p * (1 - p) / an + z2 / (4 * an**2)) / denominator
    return max(0, centre - half), min(1, centre + half)

def carrier(p):
    return 1 - (1 - p) ** 2

def panel_probability(values):
    return 1 - math.prod(1 - p for p in values)

freq = pd.read_csv(INPUT)
weight_df = pd.read_csv(WEIGHTS)
weights = dict(zip(
    weight_df[weight_df["scenario"] == "S0"]["gnomad_group"],
    weight_df[weight_df["scenario"] == "S0"]["weight"]
))

results = []

for dataset, scenario in [("exome", "S4"), ("genome", "S5")]:
    source = freq[freq["dataset"] == dataset].copy()
    source["low_an"] = source["an"] < 2500
    source["few_alt"] = source["ac"] < 20
    source["flagged"] = source["low_an"] | source["few_alt"]

    estimates = {}
    for row in source.itertuples(index=False):
        af = float(row.ac) / float(row.an)
        af_low, af_high = wilson(int(row.ac), int(row.an))
        estimates[(row.gene, row.group)] = {
            "carrier": carrier(af),
            "low": carrier(af_low),
            "high": carrier(af_high),
            "expresser": 1 - af**2,
            "expresser_low": 1 - af_high**2,
            "expresser_high": 1 - af_low**2,
            "flagged": bool(row.flagged),
        }

    quantities = {}
    for gene in GENES:
        if gene == "VKORC1" and dataset == "exome":
            quantities[gene] = None
            continue
        quantities[gene] = {
            group: estimates[(gene, group)]
            for group in weights
            if (gene, group) in estimates
        }

    for quantity in GENES + ["CYP3A5_expressers"] + list(PANELS):
        if quantity == "CYP3A5_expressers":
            genes = ["CYP3A5"]
            fields = ("expresser", "expresser_low", "expresser_high")
        elif quantity in GENES:
            genes = [quantity]
            fields = ("carrier", "low", "high")
        else:
            genes = PANELS[quantity]
            fields = ("carrier", "low", "high")

        # Exome VKORC1-dependent quantities are unavailable.
        if any(quantities.get(gene) is None for gene in genes):
            results.append({
                "quantity": quantity, "scenario": scenario,
                "national": None, "low": None, "high": None,
                "coverage": None, "weight_on_flagged_rows": None,
                "status": "N/A: VKORC1 unavailable in exome data",
            })
            continue

        common = [
            group for group in weights
            if all(group in quantities[gene] for gene in genes)
        ]
        missing = set(weights) - set(common)
        if missing:
            results.append({
                "quantity": quantity, "scenario": scenario,
                "national": None, "low": None, "high": None,
                "coverage": None, "weight_on_flagged_rows": None,
                "status": f"N/A: missing weighted groups {sorted(missing)}",
            })
            continue

        vals = []
        flagged_groups = set()
        for group in common:
            group_values = [quantities[g][group] for g in genes]
            if any(v["flagged"] for v in group_values):
                flagged_groups.add(group)

            vals.append([
                panel_probability([v[fields[0]] for v in group_values]),
                panel_probability([v[fields[1]] for v in group_values]),
                panel_probability([v[fields[2]] for v in group_values]),
            ])

        national = [
            sum(weights[g] * vals[i][j] for i, g in enumerate(common))
            for j in range(3)
        ]
        results.append({
            "quantity": quantity, "scenario": scenario,
            "national": national[0], "low": national[1], "high": national[2],
            "coverage": sum(weights[g] for g in common) * 0.949010,
            "weight_on_flagged_rows": sum(weights[g] for g in flagged_groups),
            "status": "OK",
        })

out = pd.DataFrame(results)
out.to_csv(OUTPUT, index=False)
print(out.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
print(f"\nSaved {OUTPUT}")
