from api import *


if __name__ == '__main__':
    receptor_genes = get_cpdb_receptors()
    full_ligand_set = get_ligands()
    op = process_from_web()
    op_filtered = filter_op_stmts(op.statements,
                                  full_ligand_set,
                                  receptor_genes)

    logger.info('Total OP statements %d' % (len(op_filtered)))

    # Filter statements which are not ligands/receptors from
    # OmniPath database
    op_filtered = ac.filter_direct(op_filtered)

    op_filtered = ac.filter_by_curation(op_filtered,
                                        curations=db_curations)

    # get ligands by receptor for indra_op
    op_receptor_by_ligands = get_receptor_by_ligands(receptor_genes,
                                                     full_ligand_set,
                                                     op_filtered)

    # Assemble the statements into HTML formatted report and save into a file
    html_report = \
        html_assembler(op_filtered,
                       fname=(os.path.join(OUTPUT, 'op_interactions.html')))

    nature_interactions = process_nature_paper()

    op_interactions = make_interaction_df(op_receptor_by_ligands)

    # Create a cellphone db formatted indra_op database
    op_nature = \
        [{'partner_a': i.split('_')[0],
          'partner_b': i.split('_')[1]}
         for i in set(op_interactions.interactions) | set(nature_interactions.interactions)]
    op_nature = pd.DataFrame(op_nature)
    op_nature_uniprot = []
    op_nature_interactions = []

    count = 0
    for r, c in op_nature.iterrows():
        count += 1
        if c[0] in up_hgnc and c[1] in up_hgnc:
            op_nature_uniprot.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': up_hgnc[c[0]],
                    'partner_b': up_hgnc[c[1]],
                    'source': 'INDRA'
                }
            )

            op_nature_interactions.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': c[0],
                    'partner_b': c[1],
                }
            )

    pd.DataFrame(op_nature_interactions).to_csv(os.path.join(
        HERE, os.pardir, 'output/op_nature_interactions.csv'), sep=",", index=False)

    pd.DataFrame(op_nature_uniprot).to_csv(os.path.join(
        HERE, os.pardir, 'output/op_nature_uniprot.csv'), sep=",", index=False)
