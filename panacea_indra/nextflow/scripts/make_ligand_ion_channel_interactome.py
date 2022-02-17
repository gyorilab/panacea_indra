from api import *


logger = logging.getLogger('Ion channel interactome')


full_ligand_set = get_ligands()
logger.info('Total ligands in data %d' %(len(full_ligand_set)))
ion_channels = get_ion_channels()
logger.info('Total ion channels in data %d' %(len(ion_channels)))

# Fetch omnipath database biomolecular interactions and
# process them into INDRA statements
op = process_from_web()
logger.info('Total OP statements %d' % (len(op.statements)))

# Filter statements which are not ligands/receptors from
# OmniPath database
op_filtered = filter_op_stmts(op.statements, full_ligand_set,
                              ion_channels)


# Now get INDRA DB Statements for the receptor-ligand pairs
hashes_by_gene_pair = get_hashes_by_gene_pair(indra_df, full_ligand_set,
                                              ion_channels)

# get the union of all the statement hashes
all_hashes = set.union(*hashes_by_gene_pair.values())
logger.info('Total raw statements %d' % (len(all_hashes)))

# Download the statements by hashes
stmts_by_hash = download_statements(all_hashes)
# get only the list of all the available statements
indra_db_stmts = list(stmts_by_hash.values())
logger.info('Total indra statements in data %d' %(len(indra_db_stmts)))

# Merge omnipath/INDRA statements and run assembly
indra_op_stmts = ac.run_preassembly(indra_db_stmts + op_filtered,
                                    run_refinement=False)
logger.info('Total INDRA OP statements in data %d' %(len(indra_op_stmts)))

# Filter to complex statements
indra_op_filtered = filter_to_complex_statements(indra_op_stmts,
                                                 full_ligand_set,
                                                 ion_channels)
logger.info('Statements filtered to complex %d' % (len(indra_op_filtered)))

# Filtering out the indirect INDRA statements
indra_op_filtered = ac.filter_direct(indra_op_filtered)
logger.info('Statements after filtering to direct %d' %(len(indra_op_filtered)))

# Filter by evidence
indra_op_filtered = filter_by_evidence(indra_op_filtered)
logger.info('Statements filtered to complex %d' % (len(indra_op_filtered)))






# Filter incorrect curations
indra_op_filtered = filter_incorrect_curations(indra_op_filtered)
logger.info('Total statements %d' % (len(indra_op_filtered)))

with open(os.path.join(OUTPUT, 'ion_channel_ligands_interaction_list.pkl'), 'wb') as fh:
    pickle.dump(indra_op_filtered, fh)

ion_by_ligands = get_receptor_by_ligands(ion_channels, full_ligand_set, indra_op_filtered)

# Assemble the statements into HTML formatted report and save into a file
indra_db_html_report = \
    html_assembler(indra_op_filtered,
                   fname=(os.path.join(OUTPUT, 'ligand_ion_channel_interactions.html')))
