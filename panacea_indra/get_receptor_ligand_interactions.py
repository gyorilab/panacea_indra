import tqdm
import pickle
import logging
import pandas as pd
from collections import defaultdict
from indra.util import batch_iter
from indra.databases import uniprot_client
from indra.sources import indra_db_rest
from indra.ontology.bio import bio_ontology

logger = logging.getLogger('receptor_ligand_interactions')

goa_gaf = '/Users/ben/genewalk/resources/goa_human.gaf'


def _load_goa_gaf():
    """Load the gene/GO annotations as a pandas data frame."""
    #goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
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
    #goa = goa[goa['Evidence_Code'].isin(goa_ec)]
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
    for group in tqdm.tqdm(batch_iter(hashes, 200), total=int(len(hashes)/200)):
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


if __name__ == '__main__':
    ligand_terms = ['cytokine activity', 'hormone activity',
                    'growth factor activity']
    receptor_terms = ['signaling receptor activity']

    ligand_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                     for term in ligand_terms]
    receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                       for term in receptor_terms]

    ligand_genes = get_genes_for_go_ids(ligand_go_ids)
    receptor_genes = get_genes_for_go_ids(receptor_go_ids)

    logger.info(f'Loaded {len(ligand_genes)} ligand genes')
    logger.info(f'Loaded {len(receptor_genes)} receptor genes')

    df = load_indra_df('/Users/ben/data/db_dump_df.pkl')

    hashes_by_gene_pair = get_hashes_by_gene_pair(df, ligand_genes,
                                                  receptor_genes)
    all_hashes = set.union(*hashes_by_gene_pair.values())
    stmts_by_hash = download_statements(all_hashes)
