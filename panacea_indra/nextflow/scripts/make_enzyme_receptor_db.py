import os
import logging
import pandas as pd
from make_ligand_receptor_database import get_all_enzymes


logger = logging.getLogger('Enzyme Product Interactome')

if __name__ == '__main__':
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
    logger.info('Pathway commons enzyme table fater filtering to enzymes in data: %d' % (len(pc)))

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
        if a in products and b in receptor_genes_go:
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
    stmts_to_filter = {'Complex', 'Activation', 'Inhibition'}
    boolean_series = enzyme_target_df['Interaction'].isin(stmts_to_filter)
    enzyme_target_df = enzyme_target_df[boolean_series]