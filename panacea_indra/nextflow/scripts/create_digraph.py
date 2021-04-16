import os
import sys
import pickle
import networkx
import argparse
import itertools
from collections import defaultdict
import indra.tools.assemble_corpus as ac
from indra.assemblers.html import HtmlAssembler


def create_interaction_digraph(ligand_receptors,
                               sorted_enzyme_FC,
                               fname,
                               enzyme_product_dict,
                               products_receptors):
  '''
  This function takes two dictionaries as input,
  ligand receptors and enzyme fold change and creates
  a interaction Digraph of ligands, enzymes and receptors.

  Parameters
  ----------
  celtype_stmts : Optional[list[indra.statements.Statement]]
      A list of INDRA Statements to be assembled.
  network_name : Optional[str]
      The name of the network to be assembled. Default: indra_assembled

  Attributes
  ----------
  ligands_dict : dict
      Dict of foldchange and ligands as keys and receptors as values
  enzyme dict : dict
      Dict of foldchange as keys and enzymes as values
  fname : str
      output file name
  '''

  ligand_receptors = dict(sorted(ligand_receptors.items(),
                                 reverse=True))
  G = networkx.DiGraph()

  top_lg_rc = dict(sorted(itertools.islice(ligand_receptors.items(), 30)))
  top_en = dict(itertools.islice(sorted_enzyme_FC.items(), 30))

  for FC_lg, rcs in top_lg_rc.items():
      for rc in rcs:
          G.add_node(FC_lg[1], color='green')
          G.add_edge(FC_lg[1], rc, label="{:.2f}".format(FC_lg[0]))
  for en_FC, en in top_en.items():
      for chem in enzyme_product_dict[en]:
          for rcs in products_receptors[chem]:
              G.add_node(en, color='red')
              G.add_edge(en, chem, label="{:.2f}".format(en_FC))
              G.add_edge(chem, rcs)

  G.graph.setdefault('graph', {})['rankdir'] = 'LR'
  ag = networkx.nx_agraph.to_agraph(G)
  fname = os.path.join(OUTPUT, fname + "interactions_digraph.pdf")
  ag.draw(fname, prog='dot')


def html_assembler(indra_stmts, fname):
  """Assemble INDRA statements into a HTML report"""
  html_assembler = HtmlAssembler(indra_stmts,
                                 db_rest_url='https://db.indra.bio')
  assembled_html_report = html_assembler.make_model(no_redundancy=True)
  html_assembler.save_model(fname)
  return assembled_html_report


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("--input")
  parser.add_argument("--output")
  parser.add_argument("--enzyme_product_dict")
  parser.add_argument("--products_receptors")
  parser.add_argument("--stmts_db_by_cell_type")
  parser.add_argument("--ligands_FC")
  parser.add_argument("--enzymes_FC")
  args = parser.parse_args()

  INPUT = args.input
  OUTPUT = args.output
  ENZYME_PRODUCT_DICT = args.enzyme_product_dict
  PRODUCTS_RECEPTORS = args.products_receptors
  STMT_DB_BY_CELL_TYPE = args.stmts_db_by_cell_type
  LIGANDS_FC = args.ligands_FC
  ENZYMES_FC = args.enzymes_FC


  # Di-graph of interactions between top 10 DE ligands, enzyme products
  # and receptors

  with open(ENZYME_PRODUCT_DICT, 'rb') as fh:
    enzyme_product_dict = pickle.load(fh)

  with open(PRODUCTS_RECEPTORS, 'rb') as fh:
    products_receptors = pickle.load(fh)
      
  with open(LIGANDS_FC, 'rb') as fh:
    ligands_FC = pickle.load(fh)

  with open(ENZYMES_FC, 'rb') as fh:
    enzymes_FC = pickle.load(fh)

  with open(STMT_DB_BY_CELL_TYPE, 'rb') as fh:
    stmts_db_by_cell_type = pickle.load(fh)

  with open(os.path.join(OUTPUT, 'receptors_in_data.pkl'), 'rb') as fh:
    receptors_in_data = pickle.load(fh)


  ligands_logFC = defaultdict(set)
  ligandsFC_by_receptor = defaultdict(set)
  sorted_ligands_FC = dict(sorted(ligands_FC.items(), reverse=True))
  sorted_enzyme_FC = dict(sorted(enzymes_FC.items(), reverse=True))
  lg_fc = list(sorted_ligands_FC.keys())
  lg = list(sorted_ligands_FC.values())
  filtered_stmts_by_cell_type = []
  cell_type_stmts = defaultdict(set)

  for cell_type in stmts_db_by_cell_type.keys():
    stmts = stmts_db_by_cell_type[cell_type]
    cell_type_stmts = defaultdict(set)
    for stmt in stmts:
        agent_names = {agent.name for agent in stmt.agent_list()}
        receptors = agent_names & receptors_in_data
        ligands = agent_names & set(lg)
        if len(ligands) > 0 and ligands not in receptors:
            for ligand in ligands:
                filtered_stmts_by_cell_type.append(stmt)
                lg_logFC = lg_fc[lg.index(ligand)]
                for receptor in receptors:
                    # Storing the interactions in a dictionary
                    # for each cell type and plot di-graphs
                    cell_type_stmts[(lg_logFC, ligand)].add(receptor)

                    # Storing interactions for all the cell-types
                    ligandsFC_by_receptor[(lg_logFC, ligand)].add(receptor)

    # Plot interactions for each cell type
    create_interaction_digraph(cell_type_stmts,
                               sorted_enzyme_FC,
                               os.path.join(cell_type,cell_type + "_"),
                               enzyme_product_dict,
                               products_receptors)


'''
  filtered_stmts_by_cell_type = ac.run_preassembly(filtered_stmts_by_cell_type,
                                                   run_refinement=False)
  # Assemble the statements into HTML formatted report and save into a file
  indra_op_html_report = \
      html_assembler(
          filtered_stmts_by_cell_type,
          fname=os.path.join(OUTPUT,
                             'filtered_ligand_enzyme_interactions.html'))

  # Plotting interaction Di-graph for all the cell-types
  create_interaction_digraph(ligandsFC_by_receptor, 
                             sorted_enzyme_FC,
                             '', 
                             enzyme_product_dict,
                             products_receptors)
'''