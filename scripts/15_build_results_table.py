from pathlib import Path

import pandas as pd

DATA = Path("data")
OUT = Path("results/results_table.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

genes = pd.read_csv(DATA / "carrier_by_group.csv")
panels = pd.read_csv(DATA / "panel_by_group.csv")
national = pd.read_csv(DATA / "national_estimates.csv")
sensitivity = pd.read_csv(DATA / "sensitivity_summary.csv")

QUANTITIES = [
    "ABCG2", "VKORC1", "CYP3A5", "CYP3A5_expressers", "CYP4F2", "NUDT15",
    "all_five_panel", "reduced_panel",
]
COLUMNS = [
    "quantity", "group", "estimate", "low", "high", "an", "flagged",
    "weight_on_flagged_rows", "coverage", "sensitivity_low", "sensitivity_high",
]
CANADA = "Canada (weighted)"

rows = []

# Gene rows: carrier probability for every gene, plus expressers for CYP3A5.
for _, r in genes.iterrows():
    flagged = bool(r["low_an"]) or bool(r["few_alt"])
    base = {"group": r["group"], "an": r["an"], "flagged": flagged}

    assert abs(r["carrier"] - (1 - (1 - r["af"]) ** 2)) < 1e-9, "carrier formula mismatch"
    rows.append({**base, "quantity": r["gene"], "estimate": r["carrier"],
                 "low": r["carrier_low"], "high": r["carrier_high"]})

    if r["gene"] == "CYP3A5":
        assert abs(r["expresser"] - (1 - r["af"] ** 2)) < 1e-9, "expresser formula mismatch"
        rows.append({**base, "quantity": "CYP3A5_expressers", "estimate": r["expresser"],
                     "low": r["expresser_low"], "high": r["expresser_high"]})

# Panel rows.
PANELS = {"all_five": "all_five_panel", "without_CYP3A5_VKORC1": "reduced_panel"}
for _, r in panels.iterrows():
    rows.append({"quantity": PANELS[r["panel"]], "group": r["group"],
                 "estimate": r["panel_carrier"], "low": r["panel_low"], "high": r["panel_high"],
                 "flagged": bool(r["any_low_an_or_few_alt"])})

# Canada rows (primary mapping) with the sensitivity range excluding S1.
primary = national[national["scenario"] == "primary"]
for _, r in primary.iterrows():
    s = sensitivity[sensitivity["quantity"] == r["quantity"]]
    assert len(s) == 1, f"No unique sensitivity row for {r['quantity']}"
    s = s.iloc[0]
    rows.append({"quantity": r["quantity"], "group": CANADA,
                 "estimate": r["national"], "low": r["low"], "high": r["high"],
                 "weight_on_flagged_rows": r["weight_on_flagged_rows"], "coverage": r["coverage"],
                 "sensitivity_low": s["min_excluding_S1"], "sensitivity_high": s["max_excluding_S1"]})

result = pd.DataFrame(rows).reindex(columns=COLUMNS)
result["quantity"] = pd.Categorical(result["quantity"], categories=QUANTITIES, ordered=True)
result["_canada"] = result["group"] == CANADA
result = result.sort_values(["_canada", "group", "quantity"]).drop(columns="_canada")
result["quantity"] = result["quantity"].astype(str)
result["an"] = result["an"].astype("Int64")

# Checks.
groups = sorted(genes["group"].unique())
assert len(groups) == 10
for g in groups:
    got = sorted(result.loc[result["group"] == g, "quantity"])
    assert got == sorted(QUANTITIES), f"{g}: {got}"
canada = result[result["group"] == CANADA]
assert sorted(canada["quantity"]) == sorted(QUANTITIES)
assert len(result) == 10 * 8 + 8, len(result)
assert result[["estimate", "low", "high"]].notna().all().all()
assert ((result["low"] <= result["estimate"] + 1e-12) & (result["estimate"] <= result["high"] + 1e-12)).all()
assert result[["estimate", "low", "high"]].apply(lambda c: c.between(0, 1).all()).all()
assert canada[["sensitivity_low", "sensitivity_high", "coverage", "weight_on_flagged_rows"]].notna().all().all()
assert result.loc[result["group"] != CANADA, ["sensitivity_low", "sensitivity_high"]].isna().all().all()
chk = canada.set_index("quantity")["estimate"]
for q, v in primary.set_index("quantity")["national"].items():
    assert abs(chk[q] - v) < 1e-12, q

result.to_csv(OUT, index=False, float_format="%.6f")
print(f"Saved: {OUT}")
print("Total rows:", len(result), "| groups:", len(groups), "| quantities:", len(QUANTITIES), "| Canada rows:", len(canada))
print("Flagged group rows:", int((result["flagged"] == True).sum()))
print("All checks passed.")
