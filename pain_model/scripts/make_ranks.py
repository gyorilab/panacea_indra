import os
import re
import pandas as pd


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

