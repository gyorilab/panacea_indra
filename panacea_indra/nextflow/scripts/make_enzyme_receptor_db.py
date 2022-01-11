from .api import *


logger = logging.getLogger('Enzyme Product Interactome')


def filter_by_evidence(stmts):
    filtered_hashes = set()
    readers = {'medscan', 'eidos', 'reach',
               'rlimsp', 'trips', 'sparser'}
    for stmt in stmts:
        sources = {ev.source_api for ev in stmt.evidence}
        evidence = len(stmt.evidence)

        if evidence < 2 and sources <= readers:
            continue
        elif sources == {'sparser'}:
            continue
        else:
            filtered_hashes.add(stmt.get_hash())
    return filtered_hashes


if __name__ == '__main__':
    # get receptors
    receptors_genes_go = get_cpdb_receptors() | get_ion_channels()

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
    stmts_hash = download_statements(set.union(set(enzyme_target_df.Statement_hash)), ev=10000)
    # filter to direct statements
    stmts = ac.filter_direct(stmts_hash.values())
    # filter incorrect curations
    stmts = filter_incorrect_curations(stmts)

    # filter out sparser hashes
    filtered_hashes = filter_by_evidence(stmts)

    filtered_enzyme_target_df = \
        enzyme_target_df[enzyme_target_df['Statement_hash'].isin(filtered_hashes)]

    # remove ATP statements
    filtered_enzyme_target_df = filtered_enzyme_target_df[filtered_enzyme_target_df['Product'] != 'ATP']

    filtered_stmts = download_statements(set(filtered_enzyme_target_df.Statement_hash))
    filtered_stmts = list(filtered_stmts.values())
    with open('../output/filtered_stmts.pkl', 'wb') as fh:
        pickle.dump(filtered_stmts, fh)

    # Filter to ion channels statements
    filtered_enzyme_target_ion_channels_df = \
        filtered_enzyme_target_df[filtered_enzyme_target_df['Receptor'].isin(get_ion_channels())]
    filtered_enzyme_target_ion_channels = download_statements(set(filtered_enzyme_target_ion_channels_df.Statement_hash))
    filtered_enzyme_target_ion_channels = list(filtered_enzyme_target_ion_channels.values())
    with open('../output/filtered_enzyme_target_ion_channels.pkl', 'wb') as fh:
        pickle.dump(filtered_enzyme_target_ion_channels, fh)

    # Filter to not ion channels statements
    filtered_enzyme_target_no_ion_channels_df = \
        filtered_enzyme_target_df[~filtered_enzyme_target_df['Receptor'].isin(get_ion_channels())]
    filtered_enzyme_target_no_ion_channels = download_statements(
        set(filtered_enzyme_target_no_ion_channels_df.Statement_hash))
    filtered_enzyme_target_no_ion_channels = list(filtered_enzyme_target_no_ion_channels.values())
    with open('../output/filtered_enzyme_target_no_ion_channels.pkl', 'wb') as fh:
        pickle.dump(filtered_enzyme_target_no_ion_channels, fh)

    #logger.info('Total final statements: %d' % (len(filtered_enzyme_target_df)))
    filtered_enzyme_target_df.to_csv(os.path.join(HERE, os.pardir, 'output/enzyme_product_target.csv'))

    enzyme_receptors = defaultdict(set)
    for v in filtered_enzyme_target_df.values:
        enzymes = v[0]
        receptor = v[3]
        for e in enzymes:
            enzyme_receptors[e].add(receptor)

    # Convert to a dataframe
    enzyme_receptors_df = []
    for k, v in enzyme_receptors.items():
        for receptor in v:
            enzyme_receptors_df.append(
                {
                    'Enzyme': k,
                    'Receptor': receptor
                }
        )
    enzyme_receptors_df = pd.DataFrame(enzyme_receptors_df)
    enzyme_receptors_df.to_csv(os.path.join(OUTPUT, 'enzyme_receptor.csv'), index=0)

    # cellphone db formatted database
    indra_op_df = []
    count = 0
    for v in enzyme_receptors_df.values:
        count += 1
        if v[0] in up_hgnc and v[1] in up_hgnc:
            indra_op_df.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': up_hgnc[v[0]],
                    'partner_b': up_hgnc[v[1]],
                    'source': 'INDRA_ENZYMES'
                }
            )

    enzyme_df = pd.DataFrame(indra_op_df)
    enzyme_df.to_csv(os.path.join(HERE, os.pardir, 'output/indra_op_enzyme_uniprot.csv'), sep=",", index=False)
