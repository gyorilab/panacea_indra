from api import *


logger = logging.getLogger('Enzyme Product Interactome')


if __name__ == '__main__':
    # get receptors
    receptors_genes = set(get_cpdb_receptors()) | set(get_nature_receptors()) \
                      | set(get_ion_channels())

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
    logger.info('Total enzyme-products in data: %d' % (len(products)))

    enzyme_target_df = []
    for a, b, stmt_type, hs, ec in zip(indra_df.agA_name,
                                       indra_df.agB_name,
                                       indra_df.stmt_type,
                                       indra_df.stmt_hash,
                                       indra_df.evidence_count):
        if a in products and b in receptors_genes:

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

    enzyme_target_df = pd.DataFrame(enzyme_target_df)
    stmts_to_filter = {'Complex', 'Activation'}
    boolean_series = enzyme_target_df['Interaction'].isin(stmts_to_filter)
    enzyme_target_df = enzyme_target_df[boolean_series]
    logger.info('Statements after filtering to Complex and Activation type: %d' % \
                len(set.union(set(enzyme_target_df.Statement_hash))))

    # download statements
    stmts_hash = download_statements(set.union(set(enzyme_target_df.Statement_hash)), ev=10000)
    # filter to direct statements
    stmts = ac.filter_direct(stmts_hash.values())
    logger.info('Statements after filtering to direct interactions: %d' % \
                len(stmts))
    # filter incorrect curations
    stmts = filter_incorrect_curations(stmts)

    # filter out hashes by statement evidence
    filtered_hashes = filter_by_evidence(stmts)

    filtered_enzyme_target_df = \
        enzyme_target_df[enzyme_target_df['Statement_hash'].isin(filtered_hashes)]
    logger.info('Total filtered hashes: %d' % (len(filtered_enzyme_target_df)))

    # remove ATP statements
    filtered_enzyme_target_df = filtered_enzyme_target_df[filtered_enzyme_target_df['Product'] != 'ATP']
    logger.info('Hashes after removing ATP statements: %d' % (len(filtered_enzyme_target_df)))

    filtered_stmts = download_statements(set(filtered_enzyme_target_df.Statement_hash))
    filtered_stmts = list(filtered_stmts.values())
    logger.info('Total statements after filtering: %d' % (len(filtered_stmts)))
    with open(os.path.join(OUTPUT, 'filtered_stmts.pkl'), 'wb') as fh:
        pickle.dump(filtered_stmts, fh)

    indra_db_html_report = \
        html_assembler(filtered_stmts,
                       fname=(os.path.join(OUTPUT, 'enzyme_receptor_interactions.html')))

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
    indra_df = []
    count = 0
    for v in enzyme_receptors_df.values:
        count += 1
        if v[0] in up_hgnc and v[1] in up_hgnc:
            indra_df.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': up_hgnc[v[0]],
                    'partner_b': up_hgnc[v[1]],
                    'source': 'INDRA_ENZYMES'
                }
            )

    enzyme_df = pd.DataFrame(indra_df)
    enzyme_df.to_csv(os.path.join(OUTPUT, 'indra_enzyme_uniprot.csv'), sep=",", index=False)
