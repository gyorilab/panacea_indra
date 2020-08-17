import os
import sys
import tqdm
import pickle
import logging
import datetime
import openpyxl
import pandas as pd
from collections import defaultdict
from indra.util import batch_iter
from indra.statements import Complex
from indra.databases import uniprot_client
from indra.sources import indra_db_rest
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.assemblers.html import HtmlAssembler
from indra.assemblers.cx.assembler import CxAssembler
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name


GO_ANNOTATIONS = '/Users/ben/genewalk/resources/goa_human.gaf'
INDRA_DB_PKL = '/Users/ben/data/db_dump_df.pkl'
DATA_SPREADSHEET = 'Neuroimmune gene list .xlsx'


logger = logging.getLogger('receptor_ligand_interactions')


mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}


def _load_goa_gaf():
    """Load the gene/GO annotations as a pandas data frame."""
    # goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
    #          'HGI', 'HEP', 'IBA', 'IBD'}
    goa = pd.read_csv(GO_ANNOTATIONS, sep='\t',
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


goa = _load_goa_gaf()


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def get_hashes_by_gene_pair(df, ligand_genes, receptor_genes):
    hashes_by_gene_pair = defaultdict(set)
    for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
        if a in ligand_genes and b in receptor_genes:
            hashes_by_gene_pair[(a, b)].add(hs)
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


def get_genes_for_go_ids(go_ids):
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


def mgi_to_hgnc_name(gene_list):
    """Convert given mouse gene symbols to HGNC equivalent symbols"""
    filtered_mgi = {mouse_gene_name_to_mgi[gene] for gene in gene_list
                    if gene in mouse_gene_name_to_mgi}
    hgnc_gene_set = set()
    for mgi_id in filtered_mgi:
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_gene_set.add(get_hgnc_name(hgnc_id))
    return hgnc_gene_set


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


def read_workbook(workbook):
    """ This function takes Excel workbook as an input and
    returns ligand and receptor gene list respectively.
    Input: Excel workbook with single(2 columns) or two sheets
    Condition: considers first column/sheet as ligand genes and second
    column/shet as receptor genes
    """
    wb = openpyxl.load_workbook(workbook)
    ligands = fix_dates(set([row[0].value for row in wb['logFC>0.25']][1:]))
    receptors = fix_dates(set([row[0].value
                               for row in wb['RPKM > 1.5 cfiber']][1:]))
    return ligands, receptors


def read_gene_list(infile, mode):
    gene_list = []
    try:
        with open(infile, mode) as FH:
            for eachGene in FH:
                gene_list.append(eachGene.strip("\n"))
        return gene_list

    except FileNotFoundError:
        sys.exit("Given file doesn't exist")


def filter_statements(stmts, ligands, receptors):
    for stmt in stmts:
        if isinstance(stmt, Complex):
            _filter_complex(stmt, ligands, receptors)
    return stmts


def _filter_complex(stmt, lg, rg):
    """Filter out the genes from Complex statements which
    are not present in the given ligand/receptor list"""
    stmt.members = [agent for agent in stmt.members
                    if agent.name in lg or agent.name in rg]
    return stmt


def html_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a HTML report"""
    html_assembler = HtmlAssembler(indra_stmts,
                                   db_rest_url='https://db.indra.bio')
    assembled_html_report = html_assembler.make_model()
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


def set_wd(x):
    """Set working directory to the provided path"""
    try:
        os.chdir(x)
        print("Working directory set to: "+os.getcwd())
    except FileNotFoundError:
        sys.exit("Please provide a working path.")


def read_stmts(infile, mode):
    with open(infile, mode) as FH:
        return pickle.load(FH)


def write_stmts(stmts_hash, file_name):
    with open(file_name, "wb") as FH:
        pickle.dump(stmts_hash, FH, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    # Set current working directory
    #set_wd("/Users/sbunga/PycharmProjects/INDRA/ligandReceptorInteractome/")
    raw_ligand_genes, raw_receptor_genes = read_workbook(DATA_SPREADSHEET)

    ligand_genes = mgi_to_hgnc_name(raw_ligand_genes)
    receptor_genes = mgi_to_hgnc_name(raw_receptor_genes)
    ligand_terms = ['cytokine activity', 'hormone activity',
                    'growth factor activity']
    receptor_terms = ['signaling receptor activity']

    ligand_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                     for term in ligand_terms]
    receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                       for term in receptor_terms]

    ligand_genes_go = get_genes_for_go_ids(ligand_go_ids)
    receptor_genes_go = get_genes_for_go_ids(receptor_go_ids)

    ligands_in_data = ligand_genes & ligand_genes_go
    receptors_in_data = receptor_genes & receptor_genes_go

    logger.info(f'Loaded {len(ligands_in_data)} ligand genes from data')
    logger.info(f'Loaded {len(receptors_in_data)} receptor genes from data')

    df = load_indra_df(INDRA_DB_PKL)
    hashes_by_gene_pair = get_hashes_by_gene_pair(df, ligands_in_data,
                                                  receptors_in_data)

    all_hashes = set.union(*hashes_by_gene_pair.values())
    stmts_by_hash = download_statements(all_hashes)

    # write stmts as pickle out
    # write_stmts(stmts_by_hash, "stmts_by_hash.pkl")

    # read statements pkl file
    # stmts_by_hash = read_stmts("stmts_by_hash.pkl", "rb")
    final_out = filter_statements(list(stmts_by_hash.values()),
                                  ligands_in_data, receptors_in_data)
    with open('ligand_receptors_indra_statemnts.pkl', 'wb') as fh:
        pickle.dump(final_out, fh)

    # Assemble the statements into HTML formatted report and save into a file
    assembled_html_report = html_assembler(final_out,
                                           fname="ligand_receptor_report.html")

    # Assemble the statements into Cytoscape networks and save the file into the disk
    # Optional: Please configure the indra config file in ~/.config/indra/config.ini with
    # NDEx credentials to upload the networks into the server
    cx_assembler_report, ndex_network_id = \
        cx_assembler(final_out, fname="ligand_receptor_report.cx")
    print(ndex_network_id)