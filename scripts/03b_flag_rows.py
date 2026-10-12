import numpy as np
import pandas as pd

AN_THRESHOLD = 2500   # rows with AN below this are flagged low_an
FEW_ALT_AC = 20       # rows with AC below this are flagged few_alt

path = 'data/gnomad_variant_frequencies_combined.csv'
d = pd.read_csv(path)

# make the script safe to re-run
d = d.drop(columns=[c for c in ['low_an', 'few_alt', 'not_observed'] if c in d.columns])

# AN of 0 or missing means no frequency can be calculated: store as blank
bad_an = d['an'].isna() | (d['an'] <= 0)
if bad_an.any():
    d.loc[bad_an, ['an', 'af']] = np.nan

d['low_an'] = d['an'].isna() | (d['an'] < AN_THRESHOLD)
d['few_alt'] = d['ac'] < FEW_ALT_AC
d['not_observed'] = d['ac'] == 0

d.to_csv(path, index=False)
print(d.shape)
print('AN set to blank:', int(bad_an.sum()))
print('low_an:', int(d['low_an'].sum()), '| few_alt:', int(d['few_alt'].sum()), '| not_observed:', int(d['not_observed'].sum()))
print(d[d['low_an'] | d['few_alt']][['gene', 'rsID', 'group', 'ac', 'an', 'af', 'low_an', 'few_alt']].to_string(index=False))
