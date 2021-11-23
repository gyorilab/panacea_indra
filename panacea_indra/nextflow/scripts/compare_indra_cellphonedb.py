import os
import logging
import pandas as pd

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('Compare_interactions')


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
