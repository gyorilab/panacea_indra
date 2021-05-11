import os
import sys
import tqdm
import pickle
import logging
import openpyxl
import argparse
import pandas as pd
from indra.util import batch_iter
from collections import defaultdict
from indra.statements import Complex
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from indra.assemblers.html import HtmlAssembler
from indra.sources.omnipath import process_from_web
from indra.assemblers.cx.assembler import CxAssembler
from indra_db.client.principal.curation import get_curations


logger = logging.getLogger('receptor_ligand_interactions')

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


def filter_op_stmts(op_stmts, lg, rg, nature_df):
    """ Filter out the statements which are not ligand and receptor """
    logger.info(f'Filtering {len(op_stmts)} to ligand-receptor interactions')
    filtered_stmts = [stmt for stmt in op_stmts if
                      (any(a.name in lg for a in stmt.agent_list())
                       and any(a.name in rg for a in stmt.agent_list()))]
    logger.info(f'{len(filtered_stmts)} left after filter')

    # Creating dictionary for ligands and receptors
    # from 2015 paper
    pairs = defaultdict(set)
    for r,c in nature_df.iterrows():
        pairs[(c[0])].add(c[1])

    count = 0
    not_in_paper = []

    for stmts in filtered_stmts:
        there = 0
        for stmt in stmts.agent_list():
            if stmt.name in pairs:
                for v in pairs[stmt.name]:
                    ag = [s.name for s in stmts.agent_list()]
                    if v in ag:
                        there = 1
        if there != 1:
            not_in_paper.append(stmts)
    
    return not_in_paper


def filter_incorrect_curations(stmts):
    # Filter incorrect curations
    indra_op_filtered = ac.filter_by_curation(stmts,
                                              curations=db_curations)
    return indra_op_filtered


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


def filter_out_medscan(stmts):
    logger.info('Filtering out medscan evidence on %d statements' % len(stmts))
    new_stmts = []
    for stmt in stmts:
        new_evidence = [e for e in stmt.evidence if e.source_api != 'medscan']
        if not new_evidence:
            continue
        stmt.evidence = new_evidence
        if not stmt.evidence:
            continue
        new_stmts.append(stmt)
    logger.info('%d statements after filter' % len(new_stmts))
    return new_stmts


def filter_db_only(stmts):
    new_stmts = []
    for stmt in stmts:
        sources = {ev.source_api for ev in stmt.evidence}
        if sources <= {'reach', 'sparser', 'trips', 'rlimsp', 'medscan', 'eidos'}:
            continue
        new_stmts.append(stmt)
    return new_stmts


def html_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a HTML report"""
    html_assembler = HtmlAssembler(indra_stmts,
                                   db_rest_url='https://db.indra.bio')
    assembled_html_report = html_assembler.make_model(grouping_level='statement',
                                                      no_redundancy=True)
    html_assembler.save_model(fname)
    return assembled_html_report


def cx_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a CX report"""
    cx_assembler = CxAssembler(indra_stmts)
    assembled_cx_report = cx_assembler.make_model()
    cx_assembler.save_model(fname)
    ndex_network_id = cx_assembler.upload_model(ndex_cred=None,
                                                private=True, style='default')
    return assembled_cx_report, ndex_network_id


def get_ligands_by_receptor(receptors_in_data, ligands_in_data, stmts):
    ligands_by_receptor = defaultdict(set)
    for stmt in stmts:
        agent_names = {agent.name for agent in stmt.agent_list()}
        receptors = agent_names & receptors_in_data
        ligands = agent_names & ligands_in_data
        for receptor in receptors:
            ligands_by_receptor[receptor] |= ligands
    return dict(ligands_by_receptor)


def process_df(workbook):
    wb = openpyxl.load_workbook(workbook)
    df = {
    'ligands': [row[1].value for row in wb['All.Pairs']][1:], 
    'receptors': [row[3].value for row in wb['All.Pairs']][1:]
    }
    lg_rg = pd.DataFrame(df)
    return lg_rg


db_curations = get_curations()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--output")
    #parser.add_argument("--hashes_by_gene_pair", nargs='+')
    parser.add_argument("--ligands_in_data", nargs='+')
    parser.add_argument("--enzyme_possible_drug_targets")
    parser.add_argument("--nature_paper_hashes", nargs='+')
    args = parser.parse_args()

    #HASHES_BY_GENE_PAIR = args.hashes_by_gene_pair
    LIGANDS_IN_DATA = args.ligands_in_data
    POSSIBLE_DRUG_TARGETS = args.enzyme_possible_drug_targets
    INPUT = args.input
    OUTPUT = args.output
    NATURE_HASHES = args.nature_paper_hashes



    
    stmts_by_cell_type = {}
    stmts_db_by_cell_type = {}
    possible_db_drug_targets = set()
    all_ranked_lg_df = pd.DataFrame(columns=['Interaction statement',
                                         'logFC'])


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

    with open(os.path.join(OUTPUT, 'receptors_in_data.pkl'), 'rb') as fh:
        receptors_in_data = pickle.load(fh)

    with open(POSSIBLE_DRUG_TARGETS, 'rb') as fh:
        possible_drug_targets = pickle.load(fh)

    # Get 2015 ligand receptor direct interactions dataframe from
    # the spreadsheet
    lg_rg = process_df(os.path.join(INPUT, 'ncomms8866-s3.xlsx'))


    count = 0

    op = process_from_web()

    for cell_type in IMMUNE_CELLTYPE_LIST:

        output_dir = os.path.join(OUTPUT, cell_type)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        #with open(HASHES_BY_GENE_PAIR[count], 'rb') as fh:
        #    hashes_by_gene_pair = pickle.load(fh)

        with open(LIGANDS_IN_DATA[count], 'rb') as fh:
            ligands_in_data = pickle.load(fh)

        with open(NATURE_HASHES[count], 'rb') as fh:
            hashes = pickle.load(fh)

        # get the union of all the statement hashes
        all_hashes = set.union(*hashes.values())
        # Download the statements by hashes
        stmts_by_hash = download_statements(all_hashes)
        # get only the list of all the available statemtns
        indra_db_stmts = list(stmts_by_hash.values())
        # Filtering out the indirect INDRA statements
        indra_db_stmts = ac.filter_direct(indra_db_stmts)
        # Filter statements which are not ligands/receptors from
        # OmniPath database and filter op statemeents which are not
        # in 2015 paper
        op_filtered = filter_op_stmts(op.statements, ligands_in_data.values(),
                                      receptors_in_data, lg_rg)

        # Merge omnipath/INDRA statements and run assembly
        indra_op_stmts = ac.run_preassembly(indra_db_stmts + op_filtered,
                                            run_refinement=False)

        # Filter incorrect curations
        indra_op_filtered = filter_incorrect_curations(indra_op_stmts)

        # Filter complex statements
        indra_op_filtered = filter_complex_statements(indra_op_filtered,
                                                      ligands_in_data.values(),
                                                      receptors_in_data)

        # We do this again because when removing complex members, we
        # end up with more duplicates
        indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                               run_refinement=False)
        # Filter out medscan statements
        stmts_public = filter_out_medscan(indra_op_filtered)

        # Save indra ligand receptor statements for the
        # respective cell type
        with open(os.path.join(output_dir,
                               'indra_ligand_receptor_statements.pkl'), 'wb') as fh:
            pickle.dump(stmts_public, fh)


        # Assemble the statements into HTML formatted report and save into a file
        indra_op_html_report = \
            html_assembler(
                stmts_public,
                fname=os.path.join(output_dir,
                                   'indra_ligand_receptor_report.html'))

        # Store the cell type specific indra statements in a dictionary
        stmts_by_cell_type[cell_type] = indra_op_filtered

        # Filter statements to database only and store them
        # in a separate dictionary
        stmts_db_by_cell_type[cell_type] = filter_db_only(stmts_public)

        # create a dictionary of receptors as keys and its repective
        # ligands as values
        ligands_by_receptor = get_ligands_by_receptor(receptors_in_data,
                                                      set(ligands_in_data.values()),
                                                      stmts_by_cell_type[cell_type])

        with open(cell_type+'_ligands_by_receptor.pkl', 'wb') as fh:
            pickle.dump(ligands_by_receptor, fh)

        # create a dictionary of receptors as keys and its repective
        # ligands as values from ligands and receptors with filtered database
        # statements
        ligands_by_receptor_db = get_ligands_by_receptor(receptors_in_data,
                                                         set(ligands_in_data.values()),
                                                         stmts_db_by_cell_type[cell_type])

        with open(cell_type+'_ligands_by_receptor_db.pkl', 'wb') as fh:
            pickle.dump(ligands_by_receptor_db, fh)


        # Take a union of receptors from ligands by receptor dictionary
        possible_drug_targets |= set(ligands_by_receptor.keys())
        possible_db_drug_targets |= set(ligands_by_receptor_db.keys())

        count += 1

    with open('stmts_by_cell_type.pkl', 'wb') as fh:
        pickle.dump(stmts_by_cell_type, fh)

    with open('stmts_db_by_cell_type.pkl', 'wb') as fh:
        pickle.dump(stmts_db_by_cell_type, fh)

    with open('possible_drug_targets.pkl', 'wb') as fh:
        pickle.dump(possible_drug_targets, fh)

    with open('possible_db_drug_targets.pkl', 'wb') as fh:
        pickle.dump(possible_db_drug_targets, fh)

        '''
        # Creating a dict of logFC as key and
        # ligand as its value for the respective cell type
        # from the hashes of ligands and receptors
        lg_logFC = {lg_fc[2]: lg_fc[0]
                    for lg_fc in hashes_by_gene_pair.keys()}

        logFC_stmts = defaultdict(set)
        fc_list = list(lg_logFC.keys())
        lg_list = list(lg_logFC.values())


        # Create a dict of logFC as key and statements as values
        for stmt in indra_op_filtered:
            for ag in stmt.agent_list():
                if ag.name in lg_list:
                    fc = fc_list[lg_list.index(ag.name)]
                    logFC_stmts[(fc)].add(stmt)


        # Sorting the values in descending order
        sorted_lg_stmts = dict(sorted(logFC_stmts.items(),
                                      reverse=True))

        # create a dataframe of ranked statements
        sorted_stmts_df = [
            {
                'Interaction statement': stmt,
                'logFC': fc
            }
            for fc, stmts in sorted_lg_stmts.items()
            for stmt in stmts
        ]
        ranked_df = pd.DataFrame(sorted_stmts_df)
        # write the dataframe to a TSV file
        ranked_df.to_csv(os.path.join(output_dir, cell_type + "_ligand_receptor_ranked_stmts.tsv"),
                         sep="\t", header=True, index=False)

        # Concat all the ligands receptors interaction statements
        # from all the cell types
        all_ranked_lg_df = pd.concat([all_ranked_lg_df, ranked_df],
                                     ignore_index=False)

        # unpack the set of ranked statemetns into a list
        # for assembling into a html file
        sorted_lg_stmts_list = [stmt for stmts in sorted_lg_stmts.values()
                                for stmt in stmts]

        # Assemble the statements into HTML formatted report and save into a file
        indra_op_html_report = \
            html_assembler(
                sorted_lg_stmts_list,
                fname=os.path.join(output_dir,
                                   'indra_ligand_receptor_report.html'))

        # Assemble the statements into Cytoscape networks and save the file
        # into the disk
        # Optional: Please configure the indra config file in
        # ~/.config/indra/config.ini with NDEx credentials to upload the
        # networks into the server
        
        indra_op_cx_report, ndex_network_id = \
            cx_assembler(
                stmts_by_cell_type[cell_type],
                fname=os.path.join(output_dir,
                                   'indra_ligand_receptor_report.cx'))
        

        





    with open("all_ligand_receptor_statements.pkl", 'wb') as fh:
        pickle.dump(all_ranked_lg_df, fh)

    all_ranked_lg_df.to_csv(os.path.join(OUTPUT, 'all_ligand_receptor_statements.csv'),
                            header=True,
                            index=False)
    '''




