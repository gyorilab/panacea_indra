import os
import tqdm
import pyobo
import obonet
import pickle
import logging
import openpyxl
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from indra.util import batch_iter
from collections import defaultdict
from indra.statements import Complex
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from matplotlib_venn import venn2, venn3
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.assemblers.html import HtmlAssembler
from indra.sources.omnipath import process_from_web
from indra.assemblers.cx.assembler import CxAssembler
from indra.databases import uniprot_client, hgnc_client
from indra_db.client.principal.curation import get_curations
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name

logger = logging.getLogger('receptor_ligand_interactions')

mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}

up_hgnc = {v: k for k, v in um.uniprot_gene_name.items()
           if k in um.uniprot_hgnc}

hgnc_up = {v: k for k, v in up_hgnc.items()}

db_curations = get_curations()

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, os.pardir, os.pardir, os.pardir, 'input')
OUTPUT = os.path.join(HERE, os.pardir, os.pardir, os.pardir, 'output')
INDRA_DB_PKL = os.path.join(INPUT, 'db_dump_df.pkl')
NATURE_LIGAND_RECEPTOR_SPREADSHEET = os.path.join(INPUT, 'ncomms8866-s3.xlsx')
GO_ANNOTATIONS = os.path.join(INPUT, 'goa_human.gaf')
ION_CHANNELS = os.path.join(INPUT, 'ion_channels.txt')
SURFACE_PROTEINS_WB = os.path.join(INPUT, 'Surface Proteins.xlsx')


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


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


# Load the INDRA DB DF
indra_df = load_indra_df(INDRA_DB_PKL)


def get_hashes_by_gene_pair(df, ligand_genes, receptor_genes):
    hashes_by_gene_pair = defaultdict(set)

    for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
        if a in ligand_genes and b in receptor_genes:
            hashes_by_gene_pair[(a, b)].add(hs)
    return hashes_by_gene_pair


def download_statements(hashes, ev=100):
    """Download the INDRA Statements corresponding to a set of hashes.
    """
    stmts_by_hash = {}
    for group in tqdm.tqdm(batch_iter(hashes, 200), total=int(len(hashes) / 200)):
        idbp = indra_db_rest.get_statements_by_hash(list(group),
                                                    ev_limit=ev)
        for stmt in idbp.statements:
            stmts_by_hash[stmt.get_hash()] = stmt
    return stmts_by_hash


def get_genes_for_go_ids(go_ids):
    """Return genes that are annotated with a given go ID or its children."""
    df = goa[goa['GO_ID'].isin(set(go_ids))]
    up_ids = sorted(list(set(df['DB_ID'])))
    gene_names = [uniprot_client.get_gene_name(up_id) for up_id in up_ids]
    gene_names = {g for g in gene_names if g}
    return gene_names


def filter_complex_statements(stmts, ligands, receptors):
    filtered_stmts = []
    for stmt in stmts:
        if isinstance(stmt, Complex):
            if len(stmt.members) <= 2:
                if (any(a.name in ligands for a in stmt.members)
                        and any(a.name in receptors for a in stmt.members)):
                    filtered_stmts.append(stmt)
        else:
            filtered_stmts.append(stmt)

    return filtered_stmts


def filter_op_stmts(op_stmts, lg, rg):
    """ Filter out the statements which are not ligand and receptor """
    logger.info(f'Filtering {len(op_stmts)} to ligand-receptor interactions')
    filtered_stmts = [stmt for stmt in op_stmts if
                      (any(a.name in lg for a in stmt.agent_list())
                       and any(a.name in rg for a in stmt.agent_list()))]
    logger.info(f'{len(filtered_stmts)} left after filter')
    return filtered_stmts


def html_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a HTML report"""
    html_assembler = HtmlAssembler(indra_stmts,
                                   db_rest_url='https://db.indra.bio')
    assembled_html_report = html_assembler.make_model(no_redundancy=True, grouping_level='statement')
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


def get_receptor_by_ligands(receptors_in_data, ligands_in_data, stmts):
    receptor_by_ligands = defaultdict(set)
    for stmt in stmts:
        agent_names = {agent.name for agent in stmt.agent_list()}
        receptors = agent_names & receptors_in_data
        ligands = agent_names & ligands_in_data
        for receptor in receptors:
            receptor_by_ligands[receptor] |= ligands
    return dict(receptor_by_ligands)


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


def filter_by_evidence(stmts):
    filtered_hashes = set()
    readers = {'medscan', 'eidos', 'reach',
               'rlimsp', 'trips', 'sparser', 'isi'}
    for stmt in stmts:
        sources = {ev.source_api for ev in stmt.evidence}
        evidence = len(stmt.evidence)

        if evidence < 2 and sources <= readers:
            continue
        elif sources == {'sparser'}:
            continue
        else:
            filtered_hashes.add(stmt.get_hash())
    return filtered_hashes


def filter_incorrect_curations(stmts):
    # Filter incorrect curations
    indra_op_filtered = ac.filter_by_curation(stmts,
                                              curations=db_curations)
    return indra_op_filtered


def mgi_to_hgnc_name(gene_list):
    """Convert given mouse gene symbols to HGNC equivalent symbols"""
    filtered_mgi = {mouse_gene_name_to_mgi[gene] for gene in gene_list
                    if gene in mouse_gene_name_to_mgi}
    hgnc_gene_set = set()
    for mgi_id in filtered_mgi:
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_gene_set.add(get_hgnc_name(hgnc_id))
    return hgnc_gene_set


def expand_with_child_go_terms(terms):
    all_terms = set(terms)
    for term in terms:
        child_terms = bio_ontology.get_children('GO', term)
        all_terms |= {c[1] for c in child_terms}
    return all_terms


''''
def get_receptors():
    receptor_terms = ['signaling receptor activity']
    receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                       for term in receptor_terms]
    receptor_go_ids = expand_with_child_go_terms(receptor_go_ids)
    # Filtering out the nuclear receptors from the receptor list
    receptor_go_ids = {r for r in receptor_go_ids if
                       'receptor' in bio_ontology.get_name('GO', r) or
                       'sensor' in bio_ontology.get_name('GO', r) or
                       'channel' in bio_ontology.get_name('GO', r)}
    nuclear_receptor_go_ids = expand_with_child_go_terms(['GO:0004879'])
    receptor_genes_go = get_genes_for_go_ids(receptor_go_ids) - \
        get_genes_for_go_ids(nuclear_receptor_go_ids)
    receptor_genes_go -= {'NR2C2', 'EGF'}
    # Add ION channels to the receptor list
    ion_channels = set()
    with open(ION_CHANNELS, 'r') as fh:
        for line in fh:
            ion_channels.add(line.strip())
    receptor_genes_go |= ion_channels
    return receptor_genes_go
'''


def get_go_receptors():
    receptor_terms = ['signaling receptor activity']
    receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                       for term in receptor_terms]
    receptor_go_ids = expand_with_child_go_terms(receptor_go_ids)
    # Filtering out the nuclear receptors from the receptor list
    receptor_go_ids = {r for r in receptor_go_ids if
                       'receptor' in bio_ontology.get_name('GO', r) or
                       'sensor' in bio_ontology.get_name('GO', r) or
                       'channel' in bio_ontology.get_name('GO', r)}
    nuclear_receptor_go_ids = expand_with_child_go_terms(['GO:0004879'])
    receptor_genes_go = get_genes_for_go_ids(receptor_go_ids) - \
        get_genes_for_go_ids(nuclear_receptor_go_ids)
    receptor_genes_go -= {'NR2C2', 'EGF'}
    return receptor_genes_go


def get_ligands():
    # Read and extract cell surface proteins from CSPA DB
    wb = openpyxl.load_workbook(SURFACE_PROTEINS_WB)
    surface_protein_set = set(row[4].value for row in wb['Sheet 1']
                              if row[6].value == 'yes')
    updated_surface_protein_set = {hgnc_client.get_hgnc_name(hgnc_client.get_current_hgnc_id(g))
                                   for g in surface_protein_set} - {None}
    logger.info('Got %d surface proteins from spreadsheet' %
                len(surface_protein_set))
    ligand_terms = ['cytokine activity', 'hormone activity',
                    'growth factor activity',
                    'extracellular matrix structural constituent']
    # Getting GO id's for ligands and receptors by using
    # GO terms
    ligand_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                     for term in ligand_terms]
    ligand_go_ids = expand_with_child_go_terms(ligand_go_ids)

    # Converting GO id's to gene symbols
    ligand_genes_go = get_genes_for_go_ids(ligand_go_ids)
    manual_ligands = set()
    ligand_genes_go = updated_surface_protein_set | ligand_genes_go | manual_ligands | get_cpdb_ligands()
    ligand_genes_go = {hgnc_client.get_hgnc_name(hgnc_client.get_current_hgnc_id(g))
                       for g in ligand_genes_go} - {None}
    ligand_genes_go = ligand_genes_go - (get_cpdb_receptors() | get_ion_channels())
    return ligand_genes_go


def make_interaction_df(indra_op):
    # Make INDRA_DB interactions dataframe
    interactions = []
    for receptor, ligands in indra_op.items():
        for ligand in ligands:
            interactions.append(
                {
                    'ligands': ligand,
                    'receptors': receptor,
                    'interactions': ligand + '_' + receptor
                }
            )
    interaction_df = pd.DataFrame(interactions)
    return interaction_df


def _make_nature_df(nature_interactions):
    # Make nature interactions dataframe
    nature_interactions_df = []
    for rows, cols in nature_interactions.iterrows():
        nature_interactions_df.append(
            {
                'ligands': cols[0],
                'receptors': cols[1],
                'interactions': cols[0] + '_' + cols[1]

            }
        )
    nature_interactions_df = pd.DataFrame(nature_interactions_df)
    return nature_interactions_df


def process_df(workbook):
    wb = openpyxl.load_workbook(workbook)
    df = {
        'ligands': [row[1].value for row in wb['All.Pairs']][1:],
        'receptors': [row[3].value for row in wb['All.Pairs']][1:]
    }
    lg_rg = pd.DataFrame(df)
    return lg_rg


def process_nature_paper():
    nature_xlsx = process_df(os.path.join(INPUT, 'ncomms8866-s3.xlsx'))
    nature_lg_by_rg = defaultdict(set)

    nature_dataframe = []
    count = 0

    for rows, cols in nature_xlsx.iterrows():
        nature_lg_by_rg[(cols[0])].add(cols[1])
        count += 1
        if cols[0] in up_hgnc and cols[1] in up_hgnc:
            nature_dataframe.append(
                {
                    'id_cp_interaction': 'NATURE-' + str(count),
                    'partner_a': up_hgnc[cols[0]],
                    'partner_b': up_hgnc[cols[1]],
                    'source': 'NATURE'
                }
            )

    nature_dataframe = pd.DataFrame(nature_dataframe)

    nature_dataframe.to_csv(os.path.join(OUTPUT, 'nature_uniprot.csv'),
                            sep=",", index=False)
    return _make_nature_df(nature_xlsx)


def merge_interactions(interactions, genes_file, uniprot_file):
    df = [{'partner_a': i.split('_')[0],
           'partner_b': i.split('_')[1]}
          for i in interactions]

    interactions_hgnc = pd.DataFrame(df)
    interactions_hgnc.to_csv(os.path.join(HERE, 'output', genes_file),
                             sep=",", index=False)

    cellphonedb_df = []
    count = 0
    for r, c in interactions_hgnc.iterrows():
        count += 1
        if c[0] in up_hgnc and c[1] in up_hgnc:
            cellphonedb_df.append(
                {
                    'id_cp_interaction': 'Woolf-' + str(count),
                    'partner_a': up_hgnc[c[0]],
                    'partner_b': up_hgnc[c[1]],
                    'source': 'INDRA'
                }
            )

    cellphonedb_df = pd.DataFrame(cellphonedb_df)
    cellphonedb_df.to_csv(os.path.join(HERE, 'output', uniprot_file),
                          sep=",", index=False)


def filter_to_complex_statements(stmts, ligands, receptors):
    ''' Filter statements to only complex type '''

    readers = {'medscan', 'eidos', 'reach',
               'rlimsp', 'trips', 'sparser',
               'tees', 'geneways', 'isi'}

    filtered_stmts = []
    for stmt in stmts:
        if filter_to_bel(stmt):
            filtered_stmts.append(stmt)
            continue
        if isinstance(stmt, Complex):
            if len(stmt.members) <= 2:
                if (any(a.name in ligands for a in stmt.members)
                        and any(a.name in receptors for a in stmt.members)):

                    sources = {ev.source_api for ev in stmt.evidence}
                    evidence = len(stmt.evidence)

                    if len(sources) == 1 and 'sparser' in sources:
                        continue
                    if evidence < 2 and sources <= readers:
                        continue
                    filtered_stmts.append(stmt)
    return filtered_stmts


def filter_to_bel(stmt):
    sources = {ev.source_api for ev in stmt.evidence}
    if ('bel' in sources) and (type(stmt).__name__ == 'IncreaseAmount' or
                               type(stmt).__name__ == 'Activation'):
        for ev in stmt.evidence:
            if 'bel' in ev.source_api and ev.epistemics.get('direct') == True:
                return stmt


def get_ion_channels():
    ion_channels = set()
    with open(ION_CHANNELS, 'r') as fh:
        for line in fh:
            ion_channels.add(line.strip())
    return ion_channels


def get_cpdb_receptors():
    # Load the protein data from cellphonedb
    protein_generated_cdb = \
        pd.read_csv(os.path.join(OUTPUT, 'cellphonedb_database',
                                 'cellphonedb', 'protein_generated.csv'))
    hgnc_up = {v: k for k, v in up_hgnc.items()}
    cdb_receptor = protein_generated_cdb['uniprot'][protein_generated_cdb.receptor == True]
    cdb_receptor = {hgnc_up[c] for c in cdb_receptor if c in hgnc_up}
    return cdb_receptor


def get_cpdb_ligands():
    hgnc_up = {v: k for k, v in up_hgnc.items()}
    # Load the protein data from cellphonedb
    protein_input = pd.read_csv(os.path.join('~/.cpdb/releases/v2.0.0/protein_input.csv'))
    ligands_up_in = protein_input[(protein_input['receptor'] == False)]
    ligands_up_in = set(ligands_up_in.uniprot)
    ligands_cpdb_in = {hgnc_up[c] for c in ligands_up_in if c in hgnc_up}
    return ligands_cpdb_in


def get_nature_receptors():
    nature_interactions = pd.read_excel(os.path.join(INPUT, 'ncomms8866-s3.xlsx'),
                                        sheet_name='All.Pairs')
    nature_pairs = set(nature_interactions.loc[0:, 'Pair.Name'])
    nature_receptors = {p.split('_')[1] for p in nature_pairs}
    return nature_receptors


def make_venn_plots(gene_sets: list, set_names: list, outname: str):
    font2 = {'family': 'Ariel', 'size': 8}  # use for labels
    plt.rc('font', **font2)  # sets the default font

    if len(gene_sets) != len(set_names):
        return False
    if len(gene_sets) > 3 or len(gene_sets) <= 1:
        return False

    if len(gene_sets) == 2:
        set1 = set(gene_sets[0])
        set2 = set(gene_sets[1])
        venn2([set1, set2],
              (set_names[0], set_names[1]))
        plt.title("")
        plt.savefig(os.path.join(OUTPUT, outname), dpi=300)
        plt.close()

    elif len(gene_sets) == 3:
        set1 = set(gene_sets[0])
        set2 = set(gene_sets[1])
        set3 = set(gene_sets[2])
        venn2([set1, set2, set3],
              (set_names[0], set_names[1], set_names[2]))
        plt.title("")
        plt.savefig(os.path.join(OUTPUT, outname), dpi=300)
        plt.close()


