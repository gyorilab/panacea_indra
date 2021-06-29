import os
import re
import sys
import pickle
import argparse as ag
from indra.sources import tas
from collections import defaultdict
import indra.tools.assemble_corpus as ac
from indra.statements.agent import default_ns_order
from indra.sources.omnipath import process_from_web

'''
---Small Molecule Search---

'''


if __name__ == '__main__':

    PARSER = ag.ArgumentParser()
    PARSER.add_argument('--input')
    PARSER.add_argument('--output')
    PARSER.add_argument('--targets_by_drug')
    PARSER.add_argument('--drug_bank_pkl')
    PARSER = PARSER.parse_args()

    # Parse arguments
    INPUT = PARSER.input
    OUTPUT = PARSER.output
    TARGETS_BY_DRUG = PARSER.targets_by_drug
    DRUG_BANK_PKL = PARSER.drug_bank_pkl


    # Fetch omnipath database biomolecular interactions and
    # process them into INDRA statements
    op = process_from_web()

    ### Small molecule search
    if not os.path.exists(os.path.join(INPUT,
                                       'stmts_inhibition.pkl')):
        # Process TAS statements
        tp = tas.process_from_web()

        # Read drugbank database pickle
        with open(DRUG_BANK_PKL, "rb") as fh:
            dp = pickle.load(fh)

        # Run preassembly on a list of statements
        stmts = ac.run_preassembly(tp.statements + dp, return_toplevel=False,
                                   run_refinement=False)

        # Filter the statements to a given statement type
        stmts_inhibition = ac.filter_by_type(stmts, 'Inhibition')
        

    else:
        with open(os.path.join(INPUT, 'stmts_inhibition.pkl'), 'rb') as fh:
            stmts_inhibition = pickle.load(fh)

    targets_by_drug = defaultdict(set)

    # Create a dictionary of Drugs and targets
    for stmt in stmts_inhibition:
        drug_grounding = stmt.subj.get_grounding(
            ns_order=default_ns_order + ['CHEMBL', 'PUBCHEM', 'DRUGBANK',
                                         'HMS-LINCS'])
        targets_by_drug[(stmt.subj.name, drug_grounding)].add(stmt.obj.name)

    # Write pickle outputs
    with open(TARGETS_BY_DRUG, 'wb') as fh:
        pickle.dump(targets_by_drug, fh)