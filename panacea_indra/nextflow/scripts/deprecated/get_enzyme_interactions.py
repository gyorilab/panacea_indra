import os
import re
import sys
import tqdm
import pickle
import logging
import pandas as pd
from indra.util import batch_iter
from collections import defaultdict
from indra.statements import Complex
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from indra.assemblers.html import HtmlAssembler
from indra_db.client.principal.curation import get_curations


logger = logging.getLogger('receptor_ligand_interactions')

def filter_incorrect_curations(stmts):
  # Filter incorrect curations
  indra_op_filtered = ac.filter_by_curation(stmts,
                                            curations=db_curations)
  return indra_op_filtered


def get_enzyme_product_interactions(df, de_en_df, receptors_in_data):
  hashes_by_gene_pair = defaultdict(set)
  seen_product = set()
  product_and_fc = defaultdict(set)
  for r, c in de_en_df.iterrows():
      if c[2] not in seen_product:
          product_and_fc[c[2]].add((c[0], c[3]))
      seen_product.add(c[2])

  for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
      if a in product_and_fc and b in receptors_in_data:
          enzyme_logFC = [e for v in product_and_fc[a]
                          for e in v]
          enzyme, logFC = enzyme_logFC[0], enzyme_logFC[1]
          hashes_by_gene_pair[(a, b, enzyme, logFC)].add(hs)
  return hashes_by_gene_pair


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


def html_assembler(indra_stmts, fname):
  """Assemble INDRA statements into a HTML report"""
  html_assembler = HtmlAssembler(indra_stmts,
                                 db_rest_url='https://db.indra.bio')
  assembled_html_report = html_assembler.make_model(no_redundancy=True)
  html_assembler.save_model(fname)
  return assembled_html_report


def make_interaction_df(interaction_dict):
    interaction_list = [
        {
            'Agent_A': [stmt[0].agent_list()][0][0].name,
            'Agent_B': [stmt[0].agent_list()][0][1].name,
            'Interaction type': re.match("\w+", str(stmt[0])).group(),
            'Enzyme': stmt[1],
            'logFC': fc
        }
        for fc, stmts in interaction_dict.items()
        for stmt in stmts
        if len(stmt[0].agent_list()) > 1

    ]
    df = pd.DataFrame(interaction_list)
    df = df.sort_values(by=['logFC'],
                        ascending=False)
    return df


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def filter_complex_statements(stmts, ligands, receptors):
    for stmt in stmts:
        if isinstance(stmt, Complex):
            # Statement updated by reference here
            _filter_complex(stmt, ligands, receptors)
    return stmts


def _filter_complex(stmt, lg, rg):
    """Filter out the genes from Complex statements which
    are not present in the given ligand/receptor list"""
    stmt.members = [agent for agent in stmt.members
                    if agent.name in lg or agent.name in rg]
    return stmt


db_curations = get_curations()


if __name__ == '__main__':

  INPUT = sys.argv[1]
  OUTPUT = sys.argv[2]
  DE_ENZYME_PRODUCT_LIST = sys.argv[3]
  ALL_LIGAND_RECEPTOR_STATEMENTS = sys.argv[4]


  INDRA_DB_PKL = os.path.join(INPUT, 'db_dump_df.pkl')
  # Load the INDRA DB DF
  df = load_indra_df(INDRA_DB_PKL)


  with open(DE_ENZYME_PRODUCT_LIST, 'rb') as fh:
      de_enzyme_product_list = pickle.load(fh)


  with open(os.path.join(OUTPUT, 'receptors_in_data.pkl'), 'rb') as fh:
      receptors_in_data = pickle.load(fh)


  # Get interactions for enzyme products expressed
  # in neurons
  de_enzyme_product_hash = get_enzyme_product_interactions(df, de_enzyme_product_list,
                                                           receptors_in_data)

  # Get union of all the hashes
  all_hashes = set.union(*de_enzyme_product_hash.values())

  # prdct_logFC stores enzyme product as a key and
  # its respective enzyme and logFC as value
  prdct_logFC = defaultdict(set)
  for prdct in de_enzyme_product_hash.keys():
      prdct_logFC[(prdct[0])].add((prdct[2], prdct[3]))

  # Download the statements from hashes
  stmts_by_hash = download_statements(all_hashes)
  indra_db_stmts = list(stmts_by_hash.values())

  # Filtering out the indirect INDRA statements
  indra_db_stmts = ac.filter_direct(indra_db_stmts)
  # Filtering out incorrect statemetns
  indra_op_filtered = filter_incorrect_curations(indra_db_stmts)
  # Filtering out the complex statements
  indra_op_filtered = filter_complex_statements(indra_op_filtered,
                                            set(de_enzyme_product_list['product']),
                                            receptors_in_data)
  # We do this again because when removing complex members, we
  # end up with more duplicates
  indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                         run_refinement=False)


  # Creating a dictionary of logFC and
  # its respective statement and enzyme
  logFC_stmts = defaultdict(set)
  products_receptors = defaultdict(set)
  for stmt in indra_op_filtered:
      for ag in stmt.agent_list():
          if ag.name in prdct_logFC:
              for k in prdct_logFC[ag.name]:
                  en, fc = k[0], k[1]
              logFC_stmts[(fc)].add((stmt, en))
              chem, rc = stmt.agent_list()[0].name, stmt.agent_list()[1].name
              products_receptors[(chem)].add(rc)

  # Sort the statements 
  sorted_stmts = dict(sorted(logFC_stmts.items(), reverse=True))

  # Create a list of sorted statements
  sorted_stmts_list = [stmt[0] for stmts in sorted_stmts.values()
                   for stmt in stmts]

  # Assemble the statements into HTML formatted report and save into a file
  indra_op_html_report = \
      html_assembler(
          sorted_stmts_list,
          fname=os.path.join(OUTPUT,
                             'logfc_indra_enzyme_product_neuron_report.html'))

  # Creating a dataframe of ranked enzyme products
  # and its interactions on the neuron side
  en_prdct_interaction_df = \
      make_interaction_df(sorted_stmts)

  # Write out the dataframe to a csv
  en_prdct_interaction_df.to_csv(os.path.join(OUTPUT,
                                              "ranked_de_enzyme_products_ineractions.csv"),
                                 index=True)

  # Save Products receptors object
  with open('products_receptors.pkl', 'wb') as fh:
      pickle.dump(products_receptors, fh)

  ### Ligand receptor interactions
  # Make a dataframe of ligand
  # receptor interactions
  with open(ALL_LIGAND_RECEPTOR_STATEMENTS, 'rb') as fh:
      all_ranked_lg_df = pickle.load(fh)

  ranked_lg_dict = defaultdict(set)
  for r, c in all_ranked_lg_df.iterrows():
      ranked_lg_dict[(c[1])].add((c[0], 'NA'))

  lg_rc_interaction_df = \
      make_interaction_df(ranked_lg_dict)

  # Concat enzyme product dataframe and
  # ligand receptor interaction dataframe
  full_interaction_df = pd.concat([en_prdct_interaction_df,
                                   lg_rc_interaction_df]).sort_values(by=['logFC'],
                                                                  ascending=False)
  #Write the full ligand, enzyme -> receptor to csv
  full_interaction_df.to_csv(os.path.join(OUTPUT, 'full_ranked_interaction.csv'),
                             index=True, sep=",", header=True)