import os
import re
import pandas as pd
from functools import reduce


HERE = os.path.realpath(os.path.dirname(__file__))
INPUT = os.path.join(HERE, os.pardir, 'input')
OUTPUT = os.path.join(HERE, os.pardir, 'output')

# Read all the transcriptomics files
transcript_files = {
    'motor_neuron': os.path.join(OUTPUT, 'gene_pct_motor_neuron_clusters.csv'),
    'heart': os.path.join(OUTPUT, 'gene_pct_heart_clusters.csv'),
    'cortex': os.path.join(OUTPUT, 'gene_pct_cortex_clusters.csv'),
    'drg': os.path.join(OUTPUT, 'gene_pct_drg_clusters.csv'),
    'visual_cortex': os.path.join(OUTPUT, 'gene_pct_visual_cortex_clusters.csv')
}

# Read all DRG files
proteomics_files = {
    'protein_exp': os.path.join(OUTPUT, 'protein_exp.csv'),
}

# Read phospho data
phospho_files = {
    'phospho': os.path.join(OUTPUT, 'phospho.csv')
}
#
# Check if all the files exists
all(True for k, v in transcript_files.items() if os.path.isfile(transcript_files[k]))
all(True for k, v in proteomics_files.items() if os.path.isfile(proteomics_files[k]))

# Store all the data frames into a dictionary
all_df = {}
for k, v in transcript_files.items():
    df = pd.read_csv(v)
    df.drop(df.columns[[1, 2, 3]], axis=1, inplace=True)
    columns = [k+'_'+c for c in df.columns if c != 'MOUSE_SYMBOL']
    columns.insert(0, 'MOUSE_SYMBOL')
    df.columns = columns
    all_df[k] = df

nociceptor_columns = ['drg_PEP1', 'drg_PEP2', 'drg_NP']

# Merge all the tables into one dataframe
df_merged = reduce(lambda left, right: pd.merge(left, right, on=['MOUSE_SYMBOL'],
                                                how='outer'), all_df.values()).fillna(0)
df_merged.to_csv(os.path.join(OUTPUT, 'all_gene_pct.csv'), index=False)

# Set first column as index
df_merged.set_index('MOUSE_SYMBOL', inplace=True)
new_columns = nociceptor_columns + \
              (df_merged.columns.drop(nociceptor_columns).tolist())

# reorder the columns by getting nociceptor columns to the front
df_merged = df_merged[new_columns]

# Create a rank dataframe
rank_df = pd.DataFrame(
    {
        'MOUSE_SYMBOL': df_merged.index
    }
)
rank_df[nociceptor_columns] = df_merged[nociceptor_columns].values

rank_df.set_index('MOUSE_SYMBOL', inplace=True)
# Create dummy columns
rank_df[df_merged.columns[3:]] = 0.00

# Read the proteomics file
proteomics_df = pd.read_csv(proteomics_files['protein_exp'])
proteomics_df.index = proteomics_df.MOUSE_SYMBOL
proteomics_df.drop(['Unnamed: 0', 'Description', 'MOUSE_SYMBOL'], axis=1, inplace=True)
proteomics_df = proteomics_df.loc[proteomics_df.index.drop_duplicates(keep=False), ]


# Create a rank column with 100 as base value
rank_df[['score']] = 100.00
for genes in df_merged.index:
    # get max value from 3 nociceptor clusters
    max_val = max(df_merged.loc[genes][0:3])
    # iterate over all the clusters except nociceptors
    for clusters in df_merged.columns[3:]:
        cluster_val = df_merged.loc[genes][clusters]
        # get the difference between the nociceptor and this current
        # cluster
        final = ((max_val - cluster_val)/max_val)*100
        rank_df.at[genes, clusters] = float("%.3f" % final)

        if final > 80:
            rank_df.at[genes, 'score'] = rank_df.loc[genes]['score'] + 20
        else:
            rank_df.at[genes, 'score'] = rank_df.loc[genes]['score'] - 20

        if genes in proteomics_df.index:
            for samples in proteomics_df.columns[0:3]:
                if (proteomics_df.loc[genes, 'DRG']+1)/(proteomics_df.loc[genes, samples]+1) > 3:
                    rank_df.at[genes, 'score'] = rank_df.loc[genes]['score'] + 30
                if samples == 'Heart_131c_sn.scaled':
                    if (proteomics_df.loc[genes, 'DRG'] + 1) / (proteomics_df.loc[genes, samples] + 1) > 3:
                        rank_df.at[genes, 'score'] = rank_df.loc[genes]['score'] + 60


rank_df.sort_values(by=['score'], ascending=False, inplace=True)
# re-arrange the columns
new_cols = rank_df.columns.to_list()
new_cols = new_cols[-1:] + new_cols[:-1]
rank_df = rank_df[new_cols]
rank_df.to_csv(os.path.join(OUTPUT, 'layer_2_rank_df.csv'))
