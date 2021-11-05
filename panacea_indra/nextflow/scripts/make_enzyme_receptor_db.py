import os
import tqdm
import pyobo
import obonet
import pickle
import logging
import pandas as pd
from pathlib import Path
from indra.util import batch_iter
from collections import defaultdict
from indra.sources import indra_db_rest
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.databases import uniprot_client, hgnc_client
from make_ligand_receptor_database import get_receptors

logger = logging.getLogger('Enzyme Product Interactome')

#__file__ = '/Users/sbunga/gitHub/panacea_indra/panacea_indra/nextflow/scripts/make_enzyme_receptor_db.py'
HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, os.pardir, 'input')
OUTPUT = os.path.join(HERE, os.pardir, 'output')
INDRA_DB_PKL = os.path.join(INPUT, 'db_dump_df.pkl')
GO_ANNOTATIONS = os.path.join(INPUT, 'goa_human.gaf')
ION_CHANNELS = os.path.join(INPUT, 'ion_channels.txt')


def get_all_enzymes():
    HOME = str(Path.home())
    ec_code_path = '.obo/ec-code/ec-code.obo'
    if not os.path.exists(os.path.join(HOME, ec_code_path)):
        _ = pyobo.get_id_name_mapping('ec-code')
        obo = obonet.read_obo(os.path.join(HOME, ec_code_path))
    else:
        obo = obonet.read_obo(os.path.join(HOME, ec_code_path))
    up_nodes = set()
    for node in obo.nodes:
        if node.startswith('uniprot'):
            up_nodes.add(node[8:])
    human_ups = {u for u in up_nodes if uniprot_client.is_human(u)}
    enzymes = {uniprot_client.get_gene_name(u) for u in human_ups}
    enzymes = {g for g in enzymes if not hgnc_client.is_kinase(g)}
    enzymes = {g for g in enzymes if not hgnc_client.is_phosphatase(g)}
    logger.info(f'Filtered {len(enzymes)} enzymes in total')
    return enzymes


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def filter_sparser(hashes):
    filtered_hashes = set()
    for hash, stmts in hashes.items():
        for stmt in {stmts}:
            sources = {ev.source_api for ev in stmt.evidence}
            if len(sources) == 1 and 'sparser' in sources:
                continue
            else:
                filtered_hashes.add(hash)
    return filtered_hashes


def download_statements(hashes):
    """Download the INDRA Statements corresponding to a set of hashes.
    """
    stmts_by_hash = {}
    for group in tqdm.tqdm(batch_iter(hashes, 200), total=int(len(hashes) / 200)):
        idbp = indra_db_rest.get_statements_by_hash(list(group),
                                                    ev_limit=10)
        for stmt in idbp.statements:
            stmts_by_hash[stmt.get_hash()] = stmt
    return stmts_by_hash


# Load the INDRA DB DF
indra_df = load_indra_df(INDRA_DB_PKL)


if __name__ == '__main__':
    # get receptors
    receptors_genes_go = get_receptors()
    # Enzyme product interactions
    PC_SIF_URL = ('https://www.pathwaycommons.org/archives/PC2/v12/'
                  'PathwayCommons12.Detailed.hgnc.sif.gz')
    enzymes = get_all_enzymes()
    logger.info('Total enzymes: %d' % (len(enzymes)))

    # Reading pathway commons table
    pc = pd.read_csv(PC_SIF_URL, sep='\t', header=None)
    logger.info('Length of Pathway commons enzyme table: %d' % (len(pc)))

    pc = pc[pc[1] == 'controls-production-of']
    logger.info('Pathway commons enzyme table after filtering to controls-production-of: %d' % (len(pc)))

    boolean_series = pc[0].isin(enzymes)
    pc = pc[boolean_series]
    logger.info('Pathway commons enzyme table after filtering to enzymes in data: %d' % (len(pc)))

    # create a dictionary for enzymes and its products
    # and convert the product chebi names into standard chemical
    # names
    enzyme_product = defaultdict(set)
    enzyme_target_df = []
    for r, c in pc.iterrows():
        chebi_name = bio_ontology.get_name('CHEBI', c[2])
        enzyme_product[(chebi_name)].add(c[0])

    # search for enzyme-product target
    product_targets = defaultdict(set)
    products = set(enzyme_product.keys())
    enzyme_target_df = []
    logger.info('Total enzyme-products in data: %d' % (len(products)))

    for a, b, stmt_type, hs, ec in zip(indra_df.agA_name,
                                       indra_df.agB_name,
                                       indra_df.stmt_type,
                                       indra_df.stmt_hash,
                                       indra_df.evidence_count):
        if a in products and b in receptors_genes_go:
            product_targets[(a)].add(b)
            enzyme_target_df.append(
                {
                    'Enzyme': (enzyme_product[a]),
                    'Product': a,
                    'Interaction': stmt_type,
                    'Receptor': b,
                    'Enzyme_count': len(enzyme_product[a]),
                    'Evidence_count': ec,
                    'Statement_hash': hs
                }
            )
    logger.info('Enzyme-products receptor targets: %d' % (len(product_targets)))

    enzyme_target_df = pd.DataFrame(enzyme_target_df)
    stmts_to_filter = {'Complex', 'Activation'}
    boolean_series = enzyme_target_df['Interaction'].isin(stmts_to_filter)
    enzyme_target_df = enzyme_target_df[boolean_series]

    # download statements
    stmts_hash = download_statements(set.union(set(enzyme_target_df.Statement_hash)))
    # filter out sparser hashes
    filtered_hashes = filter_sparser(stmts_hash)
    # Filter out only sparser statements
    filtered_enzyme_target_df = \
        enzyme_target_df[enzyme_target_df['Statement_hash'].isin(filtered_hashes)]
    logger.info('Total final statements: %d' % (len(filtered_enzyme_target_df)))
    filtered_enzyme_target_df.to_csv(os.path.join(HERE, os.pardir, 'output/enzyme_product_target.csv'))

    '''
    enzyme_targets = defaultdict(set)
    for k,v in enzyme_product.items():
        for p in v:
            if p in product_targets.keys():
                for i in product_targets[p]:
                    enzyme_targets[(k)].add(i)

    logger.info('Enzyme receptor targets: %d' % (len(enzyme_targets)))
    '''