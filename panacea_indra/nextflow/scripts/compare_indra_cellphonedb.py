import os
import logging
import pandas as pd
from make_ligand_receptor_database import *
from indra.databases.uniprot_client import um

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('Compare_interactions')


up_hgnc = {v: k for k, v in um.uniprot_gene_name.items()
           if k in um.uniprot_hgnc}

hgnc_up = {v:k for k,v in up_hgnc.items()}


__file__ = '/Users/sbunga/gitHub/panacea_indra/panacea_indra/nextflow/scripts/compare_indra_cellphonedb.py'
HERE = os.path.dirname(__file__)
OUTPUT = os.path.join(HERE, os.pardir, 'output/')
indra_interactions = pd.read_csv(os.path.join(OUTPUT, 'cellphonedb_database', 'indra_op_nature', 'interaction_input.csv'))
cellphone_interactions  = pd.read_csv(os.path.join(OUTPUT, 'cellphonedb_database', 'cellphonedb', 'interaction_input.csv'))

indra_partner_interactions = {v[0]+'_'+v[1] for v in indra_interactions.values}
cellphone_partner_interactions = {v[0]+'_'+v[1] for v in cellphone_interactions.values}

# Total indra interactions
logger.info('Total interactions in INDRA, OP and nature paper: %d' % (len(indra_partner_interactions)))

# Total cellphonedb interactions
logger.info('Total interactions in cellphonedb: %d' % len(cellphone_partner_interactions))

# Common interactions between indra and cellphonedb
common_interactions = indra_partner_interactions & cellphone_partner_interactions
logger.info('Common interactions b/w INDRA and cellphonedb: %d' % len(common_interactions))

# Unique interactions to INDRA
logger.info('Unique interactions to INDRA: %d' % len(indra_partner_interactions - cellphone_partner_interactions))

# Unique interactions to cellphonedb
logger.info('Unique interactions to cellphonedb: %d' % len(cellphone_partner_interactions - indra_partner_interactions))


unique_to_cdb = list(cellphone_partner_interactions - indra_partner_interactions)

# Map back the interactions to HGNC
unique_to_cdb_genes = {"_".join([hgnc_up[i.split("_")[0]], hgnc_up[i.split("_")[1]]])
                       for i in unique_to_cdb if len(i.split("_")) == 2
                       if i.split("_")[0] in hgnc_up and i.split("_")[1] in hgnc_up}



# Get all the receptors from indra interactome
rg = get_receptors()
# Get all ligands from the indra interactome
lg = get_ligands()

# Load the protein data from cellphonedb
protein_generated_cdb = \
    pd.read_csv(os.path.join(OUTPUT, 'cellphonedb_database', 'cellphonedb', 'protein_generated.csv'))

cdb_receptor = protein_generated_cdb['uniprot'][protein_generated_cdb.receptor == True]
cdb_receptor = {hgnc_up[c] for c in cdb_receptor if c in hgnc_up}


# Receptors unique to cellphonedb
unique_rg_cdb = cdb_receptor - rg
logger.info('Receptors unique to cellphonedb curation list: %d' % len(cdb_receptor - rg))

# Receptors unique to INDRA
logger.info('Receptors unique to INDRA: %d' % len(rg - cdb_receptor))

# How many receptors in INDRA are curated as Receptor in protein_generated table from cellphonedb
non_rg = protein_generated_cdb['uniprot'][protein_generated_cdb.receptor == False]
non_rg = {hgnc_up[i] for i in non_rg if i in hgnc_up}
logger.info('Total receptors in INDRA which are not annotated as receptors == False: %d' % (len(non_rg & rg)))

# Subset the protein table with the INDRA receptors which are annotated as False
#protein_generated_cdb['uniprot']
