import os
import re
import sys
import tqdm
import pickle
import logging
import openpyxl
import pandas as pd
from collections import defaultdict
from indra.util import batch_iter
from indra.databases import uniprot_client
from indra.sources import indra_db_rest
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.assemblers.html import HtmlAssembler
from indra.assemblers.cx.assembler import CxAssembler
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name

logger = logging.getLogger('receptor_ligand_interactions')

goa_gaf = '/Users/sbunga/PycharmProjects/INDRA/ligandReceptorInteractome/goa_human.gaf'
# set global variable for indradb
os.environ["INDRA_DB_REST_URL"] = "http://db.indra.bio/"


def _load_goa_gaf():
    """Load the gene/GO annotations as a pandas data frame."""
    # goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
    #          'HGI', 'HEP', 'IBA', 'IBD'}
    goa = pd.read_csv(goa_gaf, sep='\t',
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
    mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k) for k, v in um.uniprot_gene_name.items()
                              if k in um.uniprot_mgi}
    filtered_mgi = defaultdict(set)
    for genes in gene_list:
        if genes in mouse_gene_name_to_mgi.keys():
            filtered_mgi[genes].add(mouse_gene_name_to_mgi[genes])
    hgnc_gene_list = []
    for values in filtered_mgi.values():
        mgi = "".join([str(id) for id in values])
        hgnc_id = get_hgnc_from_mouse(mgi)
        hgnc_gene_list.append(str(get_hgnc_name(hgnc_id)))
    return hgnc_gene_list


def _process_sheets(wb, sheet_names):
    if len(sheet_names) == 2:
        ligand = wb[sheet_names[0]]
        receptor = wb[sheet_names[1]]
        ligand_list = [genes[0].value for genes in ligand]
        receptor_list = [genes[0].value for genes in receptor]
        return ligand_list, receptor_list
    elif len(sheet_names) == 1 and wb[sheet_names[0]].max_column == 2:
        sheet = sheet_names[0]
        ligand_list = [i[0].value for i in wb[sheet]]
        receptor_list = [i[1].value for i in wb[sheet]]
        return ligand_list, receptor_list
    else:
        return None, None


def read_workbook(workbook):
    """ This function takes Excel workbook as an input and
    returns ligand and receptor gene list respectively.
    Input: Excel workbook with single(2 columns) or two sheets
    Condition: considers first column/sheet as ligand genes and second
    column/shet as receptor genes
    """
    wb = openpyxl.load_workbook(workbook)
    sheet_names = wb.sheetnames
    ligand_list, receptor_list = _process_sheets(wb, sheet_names)
    if str(ligand_list) and str(receptor_list) == 'None':
        sys.exit("Error parsing the workbook")
    else:
        return ligand_list, receptor_list


def read_gene_list(infile, mode):
    gene_list = []
    try:
        with open(infile, mode) as FH:
            for eachGene in FH:
                gene_list.append(eachGene.strip("\n"))
        return gene_list

    except FileNotFoundError:
        sys.exit("Given file doesn't exist")


def _filter_complex(stmt, lg, rg):
    """Filter out the genes from Complex statements which
    are not present in the given ligand/receptor list"""
    stmt.members = [agent for agent in stmt.members
                    if agent.name in lg or agent.name in rg]
    return stmt


def filter_statements(ligand_list, receptor_list, statements_by_hash):
    """Maps passed in ligand and receptor gene list to the indra statements"""
    ligands = set(ligand_list)
    receptors = set(receptor_list)
    stmts_out = []
    for stmt_hash, stmt in statements_by_hash.items():
        if re.match("Complex", str(stmt)):
            stmt = _filter_complex(stmt, ligands, receptors)
        stmt_agent_names = {agent.name for agent in stmt.agent_list()}
        if stmt_agent_names & ligands and stmt_agent_names & receptors:
            stmts_out.append(stmt)
    return stmts_out


def html_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a HTML report"""
    html_assembler = HtmlAssembler(indra_stmts)
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
    set_wd("/Users/sbunga/PycharmProjects/INDRA/ligandReceptorInteractome/")
    raw_ligand_genes, raw_receptor_genes = read_workbook('neuroimmuneGeneList.xlsx')

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

    logger.info(f'Loaded {len(ligand_genes_go)} ligand genes')
    logger.info(f'Loaded {len(receptor_genes_go)} receptor genes')

    df = load_indra_df('db_dump_df.pkl')
    hashes_by_gene_pair = get_hashes_by_gene_pair(df, ligand_genes_go,
                                                  receptor_genes_go)

    all_hashes = set.union(*hashes_by_gene_pair.values())
    stmts_by_hash = download_statements(all_hashes)

    # write stmts as pickle out
    # write_stmts(stmts_by_hash, "stmts_by_hash.pkl")

    # read statements pkl file
    # stmts_by_hash = read_stmts("stmts_by_hash.pkl", "rb")

    final_out = filter_statements(ligand_genes, receptor_genes, stmts_by_hash)
    with open('ligand_receptors_indra_statemnts.pkl', 'wb') as fh:
        pickle.dump(final_out, fh)

    # Assemble the statements into HTML formatted report and save into a file
    assembled_html_report = html_assembler(final_out, fname="ligand_receptor_report.html")

    # Assemble the statements into Cytoscape networks and save the file into the disk
    # Optional: Please configure the indra config file in ~/.config/indra/config.ini with
    # NDEx credentials to upload the networks into the server
    cx_assembler_report, ndex_network_id = cx_assembler(final_out,
                                                        fname="ligand_receptor_report.cx")
    print(ndex_network_id)