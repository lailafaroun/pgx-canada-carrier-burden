from pathlib import Path
import pandas as pd

BASE = Path("data/raw/clinpgx/variantAnnotations")
OUTPUT = Path("data/evidence_counts.csv")

GENES = ["ABCG2", "CYP3A5", "CYP4F2", "NUDT15", "VKORC1"]

# Read source files as text to preserve IDs and avoid type mismatches.
study = pd.read_csv(BASE / "study_parameters.tsv", sep="\t", dtype=str)
drug = pd.read_csv(BASE / "var_drug_ann.tsv", sep="\t", dtype=str)

for df in (study, drug):
    df["Variant Annotation ID"] = (
        df["Variant Annotation ID"].fillna("").str.strip()
    )

# Standardize population labels without inferring ancestry from descriptions.
def normalize_group(value):
    value = str(value).strip()

    if not value or value.lower() == "unknown":
        return "Unknown"
    if value == "Custom, in vitro":
        return "Non-population / in vitro"
    if value == "Multiple groups" or value.startswith("Multiple groups,"):
        return "Multiple groups"

    known = {
        "European",
        "East Asian",
        "African American/Afro-Caribbean",
        "Near Eastern",
        "Central/South Asian",
        "Sub-Saharan African",
        "Latino",
        "American",
    }
    if value in known:
        return value

    return "Other / unclassified"


study["raw_group"] = (
    study["Biogeographical Groups"].fillna("").str.strip()
)
study["group"] = study["raw_group"].map(normalize_group)

# A parameter ID can have multiple rows. Keep each ID/group pair once.
study_groups = study[
    study["Variant Annotation ID"].ne("")
][["Variant Annotation ID", "group"]].drop_duplicates()

# Restrict the primary analysis to the five genes in the project.
drug["Gene"] = drug["Gene"].fillna("").str.strip()
drug["PMID"] = drug["PMID"].fillna("").str.strip()

relevant = drug[
    drug["Gene"].isin(GENES)
    & drug["Variant Annotation ID"].ne("")
    & drug["PMID"].ne("")
].copy()

# Join by annotation ID; unmatched annotations are reported in the audit.
joined = relevant.merge(
    study_groups,
    on="Variant Annotation ID",
    how="inner",
    validate="many_to_many",
)

# Counting units:
# - one PMID/group pair for a paper count
# - one annotation ID/group pair for annotation coverage
# - one PMID/gene/group combination for gene-specific paper counts
paper_group = joined[
    ["PMID", "group"]
].drop_duplicates()

annotation_group = joined[
    ["Variant Annotation ID", "group"]
].drop_duplicates()

gene_paper_group = joined[
    ["PMID", "Gene", "group"]
].drop_duplicates()

paper_counts = (
    paper_group.groupby("group")
    .size()
    .rename("unique_pmid_group_pairs")
)

annotation_counts = (
    annotation_group.groupby("group")
    .size()
    .rename("unique_annotation_group_pairs")
)

gene_counts = (
    gene_paper_group.groupby(["group", "Gene"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=GENES, fill_value=0)
    .rename(columns=lambda gene: f"{gene}_unique_pmids")
)

result = pd.concat(
    [paper_counts, annotation_counts, gene_counts],
    axis=1,
).fillna(0).reset_index()

count_columns = [
    "unique_pmid_group_pairs",
    "unique_annotation_group_pairs",
    *[f"{gene}_unique_pmids" for gene in GENES],
]
result[count_columns] = result[count_columns].astype(int)

result["share_of_pmid_group_pairs"] = (
    result["unique_pmid_group_pairs"]
    / result["unique_pmid_group_pairs"].sum()
)

result["group_type"] = result["group"].map(
    lambda group: (
        "Unknown" if group == "Unknown"
        else "Multiple groups" if group == "Multiple groups"
        else "Non-population" if group == "Non-population / in vitro"
        else "Other" if group == "Other / unclassified"
        else "Population label"
    )
)

# Put interpretable categories first, then sort by evidence count.
priority = {
    "European": 0,
    "East Asian": 1,
    "African American/Afro-Caribbean": 2,
    "Near Eastern": 3,
    "Central/South Asian": 4,
    "Sub-Saharan African": 5,
    "Latino": 6,
    "American": 7,
    "Multiple groups": 8,
    "Unknown": 9,
    "Non-population / in vitro": 10,
    "Other / unclassified": 11,
}
result["_order"] = result["group"].map(priority).fillna(99)
result = result.sort_values(
    ["_order", "unique_pmid_group_pairs"],
    ascending=[True, False],
).drop(columns="_order")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUTPUT, index=False)

# Audit
matched_ids = set(joined["Variant Annotation ID"])
relevant_ids = set(relevant["Variant Annotation ID"])

print("EVIDENCE COUNT AUDIT")
print("Drug annotation rows for five genes:", len(relevant))
print("Unique annotation IDs for five genes:", len(relevant_ids))
print("Unique annotation IDs matched:", len(matched_ids))
print("Unique annotation IDs without study parameters:",
      len(relevant_ids - matched_ids))
print("Joined rows:", len(joined))
print("Unique PMIDs in joined records:", joined["PMID"].nunique())
print("Unique PMID/group pairs:", len(paper_group))
print("Groups in output:", len(result))
print("Unknown PMID/group pairs:",
      int(result.loc[result["group"] == "Unknown",
                     "unique_pmid_group_pairs"].sum()))
print("Output:", OUTPUT)
print("\nCounts by group:")
print(result[[
    "group",
    "unique_pmid_group_pairs",
    "unique_annotation_group_pairs",
    "group_type",
]].to_string(index=False))

assert set(GENES).issubset(set(relevant["Gene"])), (
    "At least one target gene is absent from the drug annotation source."
)
assert result["unique_pmid_group_pairs"].sum() == len(paper_group)
assert result["unique_annotation_group_pairs"].sum() == len(annotation_group)

print("\nAudit checks passed.")
