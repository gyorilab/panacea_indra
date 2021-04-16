import os
import sys
import pickle
import argparse
import networkx
import itertools
from collections import defaultdict



def plot_interaction_potential(num_interactions_by_cell_type, fname):
    labels = {
        'DCs': 'Dendritic cells',
        'Dermal Macs': 'Dermal macrophages',
        'M2a': 'Reparative macrophages (2a)',
        'M2b': 'Reparative macrophages (2b)',
        'Monocytes': 'Monocytes',
        'Resident Mac': 'Resident macrophages',
        'Mast cells': 'Mast cells'
    }
    G = networkx.DiGraph()
    for cell_type, num_int in num_interactions_by_cell_type.items():
        G.add_node(cell_type, label=labels[cell_type])
        G.add_edge(cell_type, 'Neurons', label=num_int)
    ag = networkx.nx_agraph.to_agraph(G)
    ag.draw(fname, prog='dot')



def get_cell_type_stats(stmts, ligands, receptors):
    interactome = set()
    ligand_interactions = defaultdict(set)
    for stmt in stmts:
        stmt_ligands = {a.name for a in stmt.agent_list() if
                        a.name in ligands}
        stmt_receptors = {a.name for a in stmt.agent_list() if
                          a.name in receptors}
        for ligand, receptor in itertools.product(stmt_ligands,
                                                  stmt_receptors):
            interactome.add((ligand, receptor))
            ligand_interactions[ligand].add(receptor)
    return len(interactome), ligand_interactions


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--ligands_in_data", nargs='+')
    parser.add_argument("--stmts_db_by_cell_type")
    args = parser.parse_args()

    STMTS_DB_BY_CELL_TYPE = args.stmts_db_by_cell_type
    LIGANDS_IN_DATA = args.ligands_in_data
    INPUT = args.input
    OUTPUT = args.output


    num_interactions_by_cell_type = {}
    ligand_interactions_by_cell_type = {}


    with open(os.path.join(OUTPUT, 'receptors_in_data.pkl'), 'rb') as fh:
        receptors_in_data = pickle.load(fh)

    with open(STMTS_DB_BY_CELL_TYPE, 'rb') as fh:
        stmts_db_by_cell_type = pickle.load(fh)


    #IMMUNE_CELLTYPE_LIST = ['DCs',
    #                        'Dermal Macs']



    IMMUNE_CELLTYPE_LIST = ['DCs',
                            'Dermal Macs',
                            'M2a',
                            'M2b',
                            'Monocytes',
                            'Resident Mac',
                            'Mast cells'
                            ]

    count = 0 
    for cell_type in IMMUNE_CELLTYPE_LIST:

        with open(LIGANDS_IN_DATA[count], 'rb') as fh:
            ligands_in_data = pickle.load(fh)

        num_interactions_by_cell_type[cell_type], \
        ligand_interactions_by_cell_type[cell_type] = \
            get_cell_type_stats(stmts_db_by_cell_type[cell_type],
                                ligands_in_data.values(),
                                receptors_in_data)


        count += 1

    plot_interaction_potential(num_interactions_by_cell_type,
                           os.path.join(OUTPUT,
                                        'interaction_potential.pdf'))