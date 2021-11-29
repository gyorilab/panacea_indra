import os
import re
import pandas as pd
from sklearn.preprocessing import StandardScaler


HERE = os.path.realpath(os.path.dirname(__file__))
INPUT = os.path.join(HERE, os.pardir, 'input')
OUTPUT = os.path.join(HERE, os.pardir, 'output')

# Read all the transcriptomics files
files = {
    'motor_neuron': os.path.join(OUTPUT, 'gene_pct_motor_neuron_clusters.csv'),
    'heart': os.path.join(OUTPUT, 'gene_pct_heart_clusters.csv'),
    'cortex': os.path.join(OUTPUT, 'gene_pct_cortex_clusters.csv'),
    'drg': os.path.join(OUTPUT, 'gene_pct_drg_clusters.csv'),
    'visual_cortex': os.path.join(OUTPUT, 'gene_pct_visual_cortex_clusters.csv')
}

# Check if all the files exists
all(True for k, v in files.items() if os.path.isfile(files[k]))

all_df = {}
for k, v in files.items():
    df = pd.read_csv(v)
    df.drop(df.columns[[1, 2, 3]], axis=1, inplace=True)
    columns = [k+'_'+c for c in df.columns if c != 'MOUSE_SYMBOL']
    columns.insert(0, 'MOUSE_SYMBOL')
    df.columns = columns
    all_df[k] = df

df_merged = reduce(lambda left, right: pd.merge(left, right, on=['MOUSE_SYMBOL'],
                                                how='outer'), all_df.values()).fillna(0)
df_merged.to_csv(os.path.join(OUTPUT, 'all_gene_pct.csv'), index=False)