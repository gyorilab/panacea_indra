import os
import sys
import pickle
import obonet
import logging
import datetime
import openpyxl
import pandas as pd
from pathlib import Path
from indra.sources import tas
from collections import defaultdict
import indra.tools.assemble_corpus as ac
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.sources.omnipath import process_from_web
from indra.statements.agent import default_ns_order
from indra.databases import uniprot_client, hgnc_client
from indra_db.client.principal.curation import get_curations
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name


# Parse arguments
INPUT = sys.argv[1]
OUTPUT = sys.argv[2]
TARGETS_BY_DRUG = sys.argv[3]
RECEPTORS_IN_DATA = sys.argv[4]
ALL_ENZYMES = sys.argv[5]
FULL_LIGAND_SET = sys.argv[6]


GO_ANNOTATIONS = sys.argv[7]
DATA_SPREADSHEET = sys.argv[8]
DRUG_BANK_PKL = sys.argv[9]
ION_CHANNELS = sys.argv[10]
SURFACE_PROTEINS_WB = sys.argv[11]
RECEPTOR_GENES_GO = sys.argv[12]
#LIGAND_RECEPTOR_SPREADSHEET = sys.argv[12]


logger = logging.getLogger('receptor_ligand_interactions')

mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}
db_curations = get_curations()


# Functions
def fix_dates(gene_names):
    replacements = {
        datetime.datetime(2020, 3, 7, 0, 0): 'March7',
        datetime.datetime(2020, 3, 2, 0, 0): 'March2',
        datetime.datetime(2020, 3, 4, 0, 0): 'March4',
        datetime.datetime(2020, 3, 5, 0, 0): 'March5',
        datetime.datetime(2020, 3, 6, 0, 0): 'March6',
        datetime.datetime(2020, 3, 9, 0, 0): 'March9',
        datetime.datetime(2020, 3, 8, 0, 0): 'March8',
        datetime.datetime(2020, 3, 11, 0, 0): 'Mar11',
        datetime.datetime(2020, 9, 1, 0, 0): 'Sept1',
        datetime.datetime(2020, 9, 2, 0, 0): 'Sept2',
        datetime.datetime(2020, 9, 3, 0, 0): 'Sept3',
        datetime.datetime(2020, 9, 4, 0, 0): 'Sept4',
        datetime.datetime(2020, 9, 5, 0, 0): 'Sept5',
        datetime.datetime(2020, 9, 6, 0, 0): 'Sept6',
        datetime.datetime(2020, 9, 7, 0, 0): 'Sept7',
        datetime.datetime(2020, 9, 8, 0, 0): 'Sept8',
        datetime.datetime(2020, 9, 9, 0, 0): 'Sept9',
        datetime.datetime(2020, 9, 10, 0, 0): 'Sept10',
        datetime.datetime(2020, 9, 11, 0, 0): 'Sept11',
        datetime.datetime(2020, 9, 15, 0, 0): 'Sept15',
    }
    fixed_gene_names = set()
    for gene_name in gene_names:
        if isinstance(gene_name, datetime.datetime):
            fixed_gene_names.add(replacements[gene_name])
        else:
            fixed_gene_names.add(gene_name)
    return fixed_gene_names


def mgi_to_hgnc_name(gene_list):
    """Convert given mouse gene symbols to HGNC equivalent symbols"""
    filtered_mgi = {mouse_gene_name_to_mgi[gene] for gene in gene_list
                    if gene in mouse_gene_name_to_mgi}
    hgnc_gene_set = set()
    for mgi_id in filtered_mgi:
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_gene_set.add(get_hgnc_name(hgnc_id))
    return hgnc_gene_set


def get_genes_for_go_ids(go_ids, goa):
    """Return genes that are annotated with a given go ID or its children."""
    all_go_ids = set()
    for go_id in go_ids:
        children_go_ids = {ch[1] for ch in bio_ontology.get_children('GO', go_id)}
        all_go_ids.add(go_id)
        all_go_ids |= children_go_ids
    df = goa[goa['GO_ID'].isin(all_go_ids)]
    up_ids = sorted(list(set(df['DB_ID'])))
    gene_names = [uniprot_client.get_gene_name(up_id) for up_id in up_ids]
    gene_names = {g for g in gene_names if g}
    return gene_names


def filter_nuclear_receptors(receptors_go, go_term, goa):
    # Filtering out the nuclear receptors from the receptor list
    nuclear_receptors = get_genes_for_go_ids([go_term], goa)
    # Add any others that don't have the right annotation
    nuclear_receptors |= {'NR2C2'}
    filtered_receptors_go = receptors_go - nuclear_receptors
    return filtered_receptors_go


def _load_goa_gaf(GO_ANNOTATIONS_PATH):
    """Load the gene/GO annotations as a pandas data frame."""
    # goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
    #          'HGI', 'HEP', 'IBA', 'IBD'}
    goa = pd.read_csv(GO_ANNOTATIONS_PATH, sep='\t',
                      skiprows=23, dtype=str,
                      header=None,
                      names=['DB',
                             'DB_ID',
                             'DB_Symbol',
                             'Qualifier',
                             'GO_ID',
                             'DB_Reference',
                             'Evidence_Code',
                             'With_From',
                             'Aspect',
                             'DB_Object_Name',
                             'DB_Object_Synonym',
                             'DB_Object_Type',
                             'Taxon',
                             'Date',
                             'Assigned',
                             'Annotation_Extension',
                             'Gene_Product_Form_ID'])
    goa = goa.sort_values(by=['DB_ID', 'GO_ID'])
    # Filter out all "NOT" negative evidences
    goa['Qualifier'].fillna('', inplace=True)
    goa = goa[~goa['Qualifier'].str.startswith('NOT')]
    # Filter to rows with evidence code corresponding to experimental
    # evidence
    # goa = goa[goa['Evidence_Code'].isin(goa_ec)]
    return goa




def get_all_enzymes():
    HOME = str(Path.home())
    ec_code_path = '.obo/ec-code/ec-code.obo'
    if not os.path.exists(os.path.join(HOME, ec_code_path)):
        _ = pyobo.get_id_name_mapping('ec-code')
        obo = obonet.read_obo(os.path.join(HOME, ec_code_path))
    else:
        obo = obonet.read_obo(os.path.join(HOME, ec_code_path))
    up_nodes = set()
    for node in obo.nodes:
        if node.startswith('uniprot'):
            up_nodes.add(node[8:])
    human_ups = {u for u in up_nodes if uniprot_client.is_human(u)}
    enzymes = {uniprot_client.get_gene_name(u) for u in human_ups}
    enzymes = {g for g in enzymes if not hgnc_client.is_kinase(g)}
    enzymes = {g for g in enzymes if not hgnc_client.is_phosphatase(g)}
    logger.info(f'Filtered {len(enzymes)} enzymes in total')
    return enzymes


def read_workbook(workbook):
    """ This function takes Excel workbook as an input and
    returns ligand and receptor gene list respectively.
    Input: Excel workbook with single(2 columns) or two sheets
    Condition: considers first column/sheet as ligand genes and second
    column/shet as receptor genes
    """
    ligands_sheet = 'updated list of ligands '
    receptors_sheet = 'RPKM > 1.5 cfiber'
    wb = openpyxl.load_workbook(workbook)
    ligands = fix_dates(set([row[0].value for row in wb[ligands_sheet]][1:]))
    receptors = fix_dates(set([row[0].value
                               for row in wb[receptors_sheet]][1:]))
    return ligands, receptors


if __name__ == '__main__':

    GOA = _load_goa_gaf(GO_ANNOTATIONS)
    # Read and extract cell surface proteins from CSPA DB
    wb = openpyxl.load_workbook(SURFACE_PROTEINS_WB)
    surface_protein_set = set(row[4].value for row in wb['Sheet 1']
                              if row[6].value == 'yes')

    logger.info('Got %d surface proteins from spreadsheet' %
                len(surface_protein_set))
    ligand_terms = ['cytokine activity', 'hormone activity',
                    'growth factor activity']
    receptor_terms = ['signaling receptor activity']

    # Getting GO id's for ligands and receptors by using
    # GO terms
    ligand_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                     for term in ligand_terms]
    receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                       for term in receptor_terms]

    # Converting GO id's to gene symbols
    ligand_genes_go = get_genes_for_go_ids(ligand_go_ids, GOA)
    receptor_genes_go = get_genes_for_go_ids(receptor_go_ids, GOA)
    manual_ligands = {'THBS1'}


    # remove all the receptors from the surface_protein_set
    full_ligand_set = \
        (surface_protein_set - receptor_genes_go) | ligand_genes_go | \
        manual_ligands

    # Filtering out the nuclear receptors from the receptor list
    receptor_genes_go = filter_nuclear_receptors(receptor_genes_go,
                                                 'GO:0004879', GOA)

    # Add ION channels to the receptor list
    ion_channels = set()
    with open(ION_CHANNELS, 'r') as fh:
        for line in fh:
            ion_channels.add(line.strip())
    receptor_genes_go |= ion_channels



    # Fetch omnipath database biomolecular interactions and
    # process them into INDRA statements
    op = process_from_web()

    ### Small molecule search

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


    targets_by_drug = defaultdict(set)

    # Create a dictionary of Drugs and targets
    for stmt in stmts_inhibition:
        drug_grounding = stmt.subj.get_grounding(
            ns_order=default_ns_order + ['CHEMBL', 'PUBCHEM', 'DRUGBANK',
                                         'HMS-LINCS'])
        targets_by_drug[(stmt.subj.name, drug_grounding)].add(stmt.obj.name)


    # Collect lists of receptors based on GO annotations and
    # by reading the data
    # Read list of neuro immune genes from the spread sheet
    _, raw_receptor_genes = read_workbook(DATA_SPREADSHEET)
    receptor_genes = mgi_to_hgnc_name(raw_receptor_genes)
    receptors_in_data = receptor_genes & receptor_genes_go


    all_enzymes = get_all_enzymes()

    # Write pickle outputs
    with open(TARGETS_BY_DRUG, 'wb') as fh:
        pickle.dump(targets_by_drug, fh)


    with open(os.path.join(OUTPUT, RECEPTORS_IN_DATA), 'wb') as fh:
        pickle.dump(receptors_in_data, fh)


    with open(ALL_ENZYMES, 'wb') as fh:
        pickle.dump(all_enzymes, fh)

    with open(FULL_LIGAND_SET, 'wb') as fh:
        pickle.dump(full_ligand_set, fh)

    with open(RECEPTOR_GENES_GO, 'wb') as fh:
        pickle.dump(receptor_genes_go, fh)