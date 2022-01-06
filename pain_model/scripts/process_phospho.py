import os
import pandas as pd
from collections import defaultdict

# get the working path
HERE = os.getcwd()

# Read phospho data
phospho_xl = os.path.join(HERE, os.pardir,
                          'data/Primary_mouse/Proteomics/phos_PANA_16plx_Oct2020_working_forSam.xlsx')
phospho_df = pd.read_excel(phospho_xl, sheet_name='P31_16plx_Panacea_singl&comp')

# Read enriched genes
enriched_genes = pd.read_csv(os.path.join(HERE, os.pardir, 'output/enriched_genes.csv'), header=0)
# drop the first col
enriched_genes.drop('Unnamed: 0', axis=1, inplace=True)

# Rename gene symbol to mouse symbol
phospho_df.rename({'GeneSymbol': 'MOUSE_SYMBOL'}, axis=1, inplace=True)
# Filter to the enriched genes
boolean = phospho_df.MOUSE_SYMBOL.isin(enriched_genes.MOUSE_SYMBOL)
phospho_df = phospho_df[boolean]

df = dict()

for r, c in phospho_df.iterrows():
    if c[2] not in df.keys():
        df[c[2]] = {
            'Site Position': [c[4]],
            'Motif': [c[5]],
            'Brain_131_sn scaled': [c[57]],
            'Heart_131c_sn scaled': [c[58]],
            'Spinal cord_132n_sn scaled': [c[59]],
            'DRG': [(c[60]+c[61]+c[62]+c[63])/4]
        }
    else:
        df[c[2]]['Site Position'].append(c[4])
        df[c[2]]['Motif'].append(c[5])
        df[c[2]]['Brain_131_sn scaled'].append(c[57])
        df[c[2]]['Heart_131c_sn scaled'].append(c[58])
        df[c[2]]['Spinal cord_132n_sn scaled'].append(c[59])
        df[c[2]]['DRG'].append((c[60] + c[61] + c[62] + c[63]) / 4)

df = pd.DataFrame(df).transpose()

df.to_csv(os.path.join(HERE, os.pardir, 'output', 'phospho_exp.tsv'), sep='\t')