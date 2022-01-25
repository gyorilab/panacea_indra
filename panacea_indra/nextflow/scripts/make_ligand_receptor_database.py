from api import *


if __name__ == '__main__':
    receptor_genes = get_cpdb_receptors() | get_ion_channels()

    # remove all the receptors from the surface_protein_set
    full_ligand_set = get_ligands()

    if not os.path.isfile(os.path.join(OUTPUT, 'indra_op_stmts.pkl')):
        # Now get INDRA DB Statements for the receptor-ligand pairs
        hashes_by_gene_pair = get_hashes_by_gene_pair(indra_df, full_ligand_set,
                                                      receptor_genes)

        # get the union of all the statement hashes
        all_hashes = set.union(*hashes_by_gene_pair.values())
        # Download the statements by hashes
        stmts_by_hash = download_statements(all_hashes)
        # get only the list of all the available statements
        indra_db_stmts = list(stmts_by_hash.values())

        # Filtering out the indirect INDRA statements
        indra_db_stmts = ac.filter_direct(indra_db_stmts)
        logger.info('Total statements after filtering to direct ones %d' % (len(indra_db_stmts)))

        # Filter out the statements to database only
        # indra_db_stmts = filter_db_only(indra_db_stmts)
        # logger.info('Total statements after filtering to database only %d' % (len(indra_db_stmts)))

        # Fetch omnipath database biomolecular interactions and
        # process them into INDRA statements
        op = process_from_web()
        logger.info('Total OP statements %d' % (len(op.statements)))

        # Filter statements which are not ligands/receptors from
        # OmniPath database
        op_filtered = filter_op_stmts(op.statements, full_ligand_set,
                                      receptor_genes)
        op_filtered = ac.filter_direct(op_filtered)

        op_filtered = ac.filter_by_curation(op_filtered,
                                            curations=db_curations)

        # Merge omnipath/INDRA statements and run assembly
        indra_op_stmts = ac.run_preassembly(indra_db_stmts + op_filtered,
                                            run_refinement=False)
        with open(os.path.join(OUTPUT, 'indra_op_stmts.pkl'), 'wb') as fh:
            pickle.dump(indra_op_stmts, fh)
    else:
        with open(os.path.join(OUTPUT, 'indra_op_stmts.pkl'), 'rb') as fh:
            indra_op_stmts = pickle.load(fh)

    # Filter incorrect curations        
    indra_op_filtered = filter_incorrect_curations(indra_op_stmts)

    # Filter to complex statements
    indra_op_filtered = filter_to_complex_statements(indra_op_filtered,
                                                     full_ligand_set,
                                                     receptor_genes)

    # We do this again because when removing complex members, we
    # end up with more duplicates
    indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                           run_refinement=False)

    # get ligands by receptor for indra_op
    indra_op_receptor_by_ligands = get_receptor_by_ligands(receptor_genes,
                                                           full_ligand_set,
                                                           indra_op_filtered)

    # Assemble the statements into HTML formatted report and save into a file
    indra_db_html_report = \
        html_assembler(indra_op_filtered,
                       fname=(os.path.join(HERE, os.pardir, 'output', 'indra_op_interactions.html')))

    nature_interactions = process_nature_paper()
    indra_db_interactions = make_interaction_df(indra_op_receptor_by_ligands)

    # Create a cellphone db formatted indra_op database
    indra_op_nature = \
        [{'partner_a': i.split('_')[0],
          'partner_b': i.split('_')[1]}
         for i in set(indra_db_interactions.interactions) | set(nature_interactions.interactions)]
    indra_op_nature = pd.DataFrame(indra_op_nature)
    indra_op_nature_uniprot = []
    indra_op_nature_interactions = []

    count = 0
    for r, c in indra_op_nature.iterrows():
        count += 1
        if c[0] in up_hgnc and c[1] in up_hgnc:
            indra_op_nature_uniprot.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': up_hgnc[c[0]],
                    'partner_b': up_hgnc[c[1]],
                    'source': 'INDRA'
                }
            )

            indra_op_nature_interactions.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': c[0],
                    'partner_b': c[1],
                }
            )

    pd.DataFrame(indra_op_nature_interactions).to_csv(os.path.join(
        HERE, os.pardir, 'output/indra_op_nature_interactions.csv'), sep=",", index=False)

    pd.DataFrame(indra_op_nature_uniprot).to_csv(os.path.join(
        HERE, os.pardir, 'output/indra_op_nature_uniprot.csv'), sep=",", index=False)
