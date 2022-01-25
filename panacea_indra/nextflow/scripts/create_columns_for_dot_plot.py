import pandas as pd
from api import *

means_file = pd.read_csv(os.path.join(INPUT, os.pardir, 'panacea_indra', 'cellphonedb_processor',
                                      'cpdb_output', 'incision_neuro_cellphonedb_nature', 'incision_neuro', 'means.txt'),
                         sep='\t')


cluster_interactions = means_file.columns[12:]

# remove the following clusters
# NF1, NF2, NF3, p_cLTMR, keratinocytes
remove = ['NF1', 'NF2', 'NF3', 'p_cLTMR2', 'Keratinocytes']
clusters_removed = [c for c in cluster_interactions
                    if c.split('|')[0] not in remove and c.split('|')[1] not in remove]

neurons = ["NP", "PEP1", "PEP2", "cLTMR1", "SST"]

# Dont include neuron | neuron
clusters_removed = set(clusters_removed) - {c for c in clusters_removed
                                            if(c.split('|')[0]  in neurons and c.split('|')[1] in neurons)}

# set direction to Immune cells to left and neurons to right
clusters_direction_set = {c for c in clusters_removed
                          if(c.split('|')[0] not in neurons) and
                          (c.split('|')[1] in neurons)}

df = pd.DataFrame({'clusters':list(clusters_direction_set)})
df.to_csv(os.path.join(INPUT, os.pardir, 'panacea_indra', 'cellphonedb_processor',
                       'cpdb_output', 'incision_neuro_cellphonedb_nature', 'incision_neuro', 'columns.txt'),
          index=False, header=False)
