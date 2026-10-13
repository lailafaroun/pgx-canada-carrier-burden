import math
from pathlib import Path

import pandas as pd


INPUT = Path("data/gnomad_variant_frequencies_combined.csv")
OUTPUT = Path("data/carrier_by_group.csv")


def wilson_interval(ac, an, z=1.96):
    """Return the Wilson confidence interval for allele frequency."""
    if an <= 0:
        raise ValueError("AN must be greater than zero.")
    if ac < 0 or ac > an:
        raise ValueError("AC must be between zero and AN.")

    p = ac / an
    z2 = z ** 2
    denominator = 1 + z2 / an

    centre = (p + z2 / (2 * an)) / denominator
    half_width = (
        z * math.sqrt(p * (1 - p) / an + z2 / (4 * an ** 2))
        / denominator
    )

    return max(0.0, centre - half_width), min(
        1.0, centre + half_width
    )


def carrier_probability(p):
    """Probability of at least one copy under Hardy-Weinberg equilibrium."""
    if not 0 <= p <= 1:
        raise ValueError("Allele frequency must be between 0 and 1.")
    return 1 - (1 - p) ** 2


def expresser_probability(af):
    """Estimated share with at least one non-*3 CYP3A5 allele."""
    if not 0 <= af <= 1:
        raise ValueError("Allele frequency must be between 0 and 1.")
    return 1 - af ** 2


def run_tests():
    low, high = wilson_interval(0, 1000)
    assert low == 0
    assert abs(high - 0.0038) < 0.0002

    low, high = wilson_interval(500, 1000)
    assert abs(low - 0.469) < 0.001
    assert abs(high - 0.531) < 0.001

    assert abs(carrier_probability(0.1) - 0.19) < 1e-9

    af_low, af_high = wilson_interval(100, 1000)
    exp_estimate = expresser_probability(0.1)
    exp_low = expresser_probability(af_high)
    exp_high = expresser_probability(af_low)

    assert abs(af_low - 0.0829) < 0.001
    assert abs(af_high - 0.1202) < 0.001
    assert abs(exp_estimate - 0.9900) < 0.0001
    assert abs(exp_low - 0.9856) < 0.0002
    assert abs(exp_high - 0.9931) < 0.0002
    assert exp_low <= exp_estimate <= exp_high

    print("All function tests passed.")


def main():
    run_tests()

    df = pd.read_csv(INPUT)

    required = {"gene", "ac", "an", "af"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df[["ac", "an", "af"]].isna().any().any():
        raise ValueError("Missing AC, AN, or AF values need investigation.")

    estimates = []

    for row in df.itertuples(index=False):
        ac = int(row.ac)
        an = int(row.an)
        af = float(row.af)

        if an <= 0 or not 0 <= ac <= an:
            raise ValueError(f"Invalid AC/AN: {ac}/{an}")
        if not 0 <= af <= 1:
            raise ValueError(f"Invalid AF: {af}")
        if abs(af - ac / an) > 1e-10:
            raise ValueError(
                f"AF does not match AC/AN for {row.gene}, {row.group}"
            )

        af_low, af_high = wilson_interval(ac, an)

        carrier = carrier_probability(af)
        carrier_low = carrier_probability(af_low)
        carrier_high = carrier_probability(af_high)

        expresser = expresser_low = expresser_high = float("nan")

        if row.gene == "CYP3A5":
            expresser = expresser_probability(af)
            # This function decreases as AF increases, so reverse the bounds.
            expresser_low = expresser_probability(af_high)
            expresser_high = expresser_probability(af_low)

        estimates.append({
            "af_low": af_low,
            "af_high": af_high,
            "carrier": carrier,
            "carrier_low": carrier_low,
            "carrier_high": carrier_high,
            "expresser": expresser,
            "expresser_low": expresser_low,
            "expresser_high": expresser_high,
        })

    estimate_df = pd.DataFrame(estimates, index=df.index)
    df = pd.concat([df, estimate_df], axis=1)

    df.to_csv(OUTPUT, index=False)

    print(f"Processed {len(df)} rows.")
    print(f"Saved results to {OUTPUT}")
    print("\nNUDT15 Amish check:")
    print(df[(df["gene"] == "NUDT15") & (df["group"] == "Amish")][
        ["ac", "an", "af", "af_low", "af_high",
         "carrier", "carrier_low", "carrier_high"]
    ].to_string(index=False))

    print("\nCYP3A5 expresser examples:")
    print(df[df["gene"] == "CYP3A5"][
        ["group", "af", "expresser", "expresser_low", "expresser_high"]
    ].head().to_string(index=False))


if __name__ == "__main__":
    main()
