import os
import logging
import pandas as pd
from collections import defaultdict
from make_ligand_receptor_database import *
from indra.databases.uniprot_client import um

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Compare_interactions')

up_hgnc = {v: k for k, v in um.uniprot_gene_name.items()
           if k in um.uniprot_hgnc}

hgnc_up = {v: k for k, v in up_hgnc.items()}

__file__ = '/Users/sbunga/gitHub/panacea_indra/panacea_indra/nextflow/scripts/compare_indra_cellphonedb.py'
HERE = os.path.dirname(__file__)
OUTPUT = os.path.join(HERE, os.pardir, 'output/')
indra_interactions = pd.read_csv(
    os.path.join(OUTPUT, 'cellphonedb_database', 'indra_op_nature', 'interaction_input.csv'))
cellphone_interactions = pd.read_csv(
    os.path.join(OUTPUT, 'cellphonedb_database', 'cellphonedb', 'interaction_input.csv'))

indra_partner_interactions = {v[0] + '_' + v[1] for v in indra_interactions.values}
cellphone_partner_interactions = {v[0] + '_' + v[1] for v in cellphone_interactions.values}

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
pd.DataFrame({'Receptors': list(rg)}).to_csv(os.path.join(OUTPUT, 'go_receptors.csv'),
                                             index=False, header=True)
# Get all ligands from the indra interactome
lg = get_ligands()

# Get cellphonedb receptors
cdb_receptor = get_cdb_receptors()
protein_generated_cdb = \
    pd.read_csv(os.path.join(OUTPUT, 'cellphonedb_database', 'cellphonedb', 'protein_generated.csv'))
pd.DataFrame({'Receptors': list(cdb_receptor)}).to_csv(os.path.join(OUTPUT, 'cdb_receptors.csv'),
                                                       index=False, header=True)

# Receptors unique to cellphonedb
unique_rg_cdb = cdb_receptor - rg
logger.info('Receptors unique to cellphonedb curation list: %d' % len(cdb_receptor - rg))
# save the list
pd.DataFrame({'Receptors': list(unique_rg_cdb)}).to_csv(os.path.join(OUTPUT, 'unique_cdb_receptors.csv'),
                                                        index=False, header=True)
# Receptors unique to INDRA
logger.info('Receptors unique to INDRA: %d' % len(rg - cdb_receptor))
pd.DataFrame({'Receptors': list(rg - cdb_receptor)}).to_csv(os.path.join(OUTPUT, 'go_unique_receptors.csv'),
                                                            index=False, header=True)

# How many receptors in INDRA are curated as Receptor in protein_generated table from cellphonedb
non_rg = protein_generated_cdb['uniprot'][protein_generated_cdb.receptor == False]
non_rg = {hgnc_up[i] for i in non_rg if i in hgnc_up}
logger.info('Total receptors in INDRA which are NOT curated as receptors by cellphonedb: %d' % (len(non_rg & rg)))

logger.info('Total receptors in INDRA which are curated as receptors by cellphonedb: %d' % (len(cdb_receptor & rg)))

# save receptors from INDRA which are curated as False by cellphonedb
pd.DataFrame({'Receptors': list(non_rg & rg)}).to_csv(os.path.join(OUTPUT, 'indra_not_receptors_by_cdb.csv'),
                                                      index=False, header=True)

# get ion channels
ion_channels = get_ion_channels()
indra_non_rg = set(non_rg & rg)
logger.info('Total ion channels in INDRA receptors and not in cellphonedb: %d' % (len(ion_channels & indra_non_rg)))

# Now lets check to which ontology these remaining receptors are classified into
receptor_terms = ['signaling receptor activity']
receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                   for term in receptor_terms]
receptor_go_ids = expand_with_child_go_terms(receptor_go_ids)

# Lets create a dictionary of ontology terms as keys and
# its genes as values
receptors_ontology = defaultdict(set)
for r in receptor_go_ids:
    if 'receptor' in bio_ontology.get_name('GO', r):
        g = get_genes_for_go_ids([r])
        receptors_ontology[('receptors')].update(g)
    elif 'sensor' in bio_ontology.get_name('GO', r):
        g = get_genes_for_go_ids([r])
        receptors_ontology[('sensor')].update(g)
    elif 'channel' in bio_ontology.get_name('GO', r):
        g = get_genes_for_go_ids([r])
        receptors_ontology[('channel')].update(g)

# check to which ontology these missing receptors belong to
len((indra_non_rg - ion_channels) & receptors_ontology[('receptors')])

# Subset the protein table with the INDRA receptors which are annotated as False
