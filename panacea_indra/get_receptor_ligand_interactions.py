import pickle
import logging
import pandas as pd
from indra.databases import uniprot_client

logger = logging.getLogger('receptor_ligand_interactions')

goa_gaf = '/Users/ben/genewalk/resources/goa_human.gaf'


def _load_goa_gaf():
    """Load the gene/GO annotations as a pandas data frame."""
    goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
              'HGI', 'HEP', 'IBA', 'IBD'}
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
    goa = goa[goa['Evidence_Code'].isin(goa_ec)]
    return goa


goa = _load_goa_gaf()


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def get_genes_for_go_id(go_id):
    """Return genes that are annotated with a given go ID."""
    df = goa[goa['GO_ID'] == go_id]
    up_ids = sorted(list(set(df['DB_ID'])))
    gene_names = [uniprot_client.get_gene_name(up_id) for up_id in up_ids]
    gene_names = [g for g in gene_names if g]
    return gene_names
