import os
import pickle
import pandas as pd

wd = __file__

'''
INDRA SIF location on S3:
s3://bigmech/indra-db/dumps/2021-01-26/sif.pkl
'''
INDRA_SIF = os.path.join(os.pardir, 'input', 'indra_phospho_df.pkl')
with open(INDRA_SIF, 'rb') as fh:
    SIF = pickle.load(fh)
SIF.to_csv(os.path.join(os.pardir, 'output', 'indra_phospho_df.csv'), 
           index=0)