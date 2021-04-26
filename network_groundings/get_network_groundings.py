import os
import gilda
import pickle
import logging
import pandas as pd
from indra.sources import tas
from collections import defaultdict
from indra.statements import Complex
import indra.tools.assemble_corpus as ac
from indra.sources import indra_db_rest as idr
from indra.assemblers.html import HtmlAssembler
from indra.statements.agent import default_ns_order
from indra_db.client.principal.curation import get_curations


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


def filter_complex_statements(stmts, subj, obj):
    for stmt in stmts:
        if isinstance(stmt, Complex):
            # Statement updated by reference here
            _filter_complex(stmt, subj, obj)
    return stmts


def _filter_complex(stmt, subj, obj):
    """Filter out the genes from Complex statements which
    are not present in the given ligand/receptor list"""
    stmt.members = [agent for agent in stmt.members
                    if agent.name in subj or agent.name in obj]
    return stmt

# Filter incorrect curations
db_curations = get_curations()


HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, 'input')
OUTPUT = os.path.join(HERE, 'output')
INDRA_DB_PKL = os.path.join(INPUT, 'db_dump_df.pkl')
#DRUG_BANK_PKL = os.path.join(INPUT, 'drugbank_5.1.pkl')

files = {
    'FCN159': os.path.join(INPUT, 'FCN159', 'FCN159_all relationships.xls'),
    'P22077': os.path.join(INPUT, 'P22077', 'p22077.xls'),
    'SP2509': os.path.join(INPUT, 'SP2509', 'SP2509_all relationships.xls')
}

# Load the INDRA DB DF
with open(INDRA_DB_PKL, 'rb') as fh:
    indra_db = pickle.load(fh)

    
for file in files:
    normalized_df = []
    interaction_network = defaultdict(set)
    df = pd.read_excel(files[file], skiprows=1)
    stmts = []
    subj_set = set()
    obj_set = set()
    for r, c in df.iterrows():
        # Grounding subject and object names 
        # using Gilda
        subj, obj = c[0], c[2]
        subj_set.add(subj)
        obj_set.add(obj)
        gilda_subj = gilda.ground(subj)
        gilda_subj = gilda_subj[0].term.entry_name if gilda_subj else 'NA'
        
        gilda_obj = gilda.ground(obj)
        gilda_obj = gilda_obj[0].term.entry_name if gilda_obj else 'NA'
        
        normalized_df.append({
            'Subject' : subj,
            'Normalized subject' : gilda_subj,
            'Object' : obj,
            'Normalized object' : gilda_obj
        })
        # Downloading statements using INDRA REST API
        idrp = idr.get_statements(subject=gilda_subj, object=gilda_obj)
        stmts = stmts+idrp.statements
        
    # Filtering out the indirect INDRA statements
    #indra_stmts = ac.filter_direct(stmts)
    indra_stmts = ac.run_preassembly(stmts,
                                     run_refinement=False)
    indra_filtered = ac.filter_by_curation(indra_stmts,
                                           curations=db_curations)
    

    indra_op_filtered = filter_complex_statements(indra_filtered,
                                                  subj_set,
                                                  obj_set)
    
    indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                     run_refinement=False)
    
    html_assembler(indra_op_filtered, os.path.join(OUTPUT, file+'_indra_report.html'))
        
    normalized_df = pd.DataFrame(normalized_df)
    normalized_df.to_csv(os.path.join(INPUT, file, file+'_normalized_names.csv'))