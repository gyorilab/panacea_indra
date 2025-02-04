import os
import tqdm
import pickle
import pandas as pd
from indra.util import batch_iter
from indra.sources import indra_db_rest
from indra.statements import DecreaseAmount, IncreaseAmount
from indra.tools.assemble_corpus import filter_human_only, filter_genes_only, filter_transcription_factor

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


def filter_db_only(stmts):
    new_stmts = []
    for stmt in stmts:
        sources = {ev.source_api for ev in stmt.evidence}
        if sources <= {'reach', 'sparser', 'trips', 'rlimsp', 'medscan', 'eidos'}:
            continue
        new_stmts.append(stmt)
    return new_stmts


def make_dataframe(stmts):
    tf_df = []
    for stmt in tqdm.tqdm(stmts):
        agA = stmt.subj.name
        agB = stmt.obj.name

        tf_df.append({
            'agA' : agA,
            'agB' : agB,
            'stmt_type' : type(stmt).__name__
        })
    return pd.DataFrame(tf_df)


wd = __file__

INDRA_SIF = os.path.join(os.pardir, 'input', 'sif.pkl')
with open(INDRA_SIF, 'rb') as fh:
    SIF = pickle.load(fh)
    
n_stmt_type = list(SIF.columns).index('stmt_type')
n_stmt_hash = list(SIF.columns).index('stmt_hash')
hash_set = set()
for r,c in SIF.iterrows():
    if c[n_stmt_type] == 'IncreaseAmount' or c[n_stmt_type] == 'DecreaseAmount':
        hash_set.add(c[n_stmt_hash])

#stmts = download_statements(hash_set)
indra_stmts = list(stmts.values())
with open('../output/all_stmts.pkl', 'wb') as fh:
    pickle.dump(indra_stmts, fh)
    
indra_stmts = filter_human_only(indra_stmts)
indra_stmts = filter_genes_only(indra_stmts)
indra_stmts = filter_transcription_factor(indra_stmts)
indra_stmts_db_only = filter_db_only(indra_stmts)


indra_stmts_df = make_dataframe(indra_stmts)
indra_stmts_df.to_csv('../output/indra_all_tf.csv')


indra_stmts_db_only_df = make_dataframe(indra_stmts_db_only)
indra_stmts_db_only_df.to_csv('../output/indra_db_only_tf.csv')