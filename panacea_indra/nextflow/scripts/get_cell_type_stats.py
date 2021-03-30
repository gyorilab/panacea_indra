		# Getting cell type stats
        num_interactions_by_cell_type[cell_type], \
        ligand_interactions_by_cell_type[cell_type] = \
            get_cell_type_stats(stmts_db_by_cell_type[cell_type],
                                ligands_in_data,
                                receptors_in_data)