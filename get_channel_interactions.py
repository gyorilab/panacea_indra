import os
import json
import numpy
import boto3
import pickle
import pandas
from indra.statements import stmts_to_json
from indra.assemblers.tsv import TsvAssembler
from indra.assemblers.html import HtmlAssembler
from indra.tools import assemble_corpus as ac
from indra.belief import BeliefEngine
from indra.sources.indra_db_rest import get_statements


def get_channel_statements(channels, ev_limit=100):
    """Get all statements from the database for a list of gene symbols."""
    all_statements = {}
    for channel in channels:
        idbp = get_statements(agents=[channel], ev_limit=ev_limit,
                              best_first=False)
        ev_counts = {int(k): v for k, v in idbp.get_ev_counts().items()}
        source_counts = idbp.get_source_counts()
        stmts = idbp.statements
        stmts = filter_out_medscan(stmts, source_counts)
        stmts = filter_out_other_channels(channel, stmts,
                                          all_channel_gene_names)
        all_statements[channel] = (stmts, ev_counts, source_counts)
    return all_statements


def non_medscan_evidence(stmt, source_counts):
    counts = source_counts.get(stmt.get_hash())
    ev_count = sum(c for k, c in counts.items() if k != 'medscan')
    return ev_count


def filter_out_medscan(stmts, source_counts):
    new_stmts = []
    for stmt in stmts:
        new_evidence = []
        for ev in stmt.evidence:
            if ev.source_api == 'medscan':
                continue
            new_evidence.append(ev)
        if not non_medscan_evidence(stmt, source_counts):
            continue
        new_stmts.append(stmt)
    return new_stmts


def filter_out_other_channels(channel, stmts, channels):
    new_stmts = []
    for stmt in stmts:
        agents = [a for a in stmt.agent_list() if a is not None]
        if len(agents) < 2:
            continue
        uniq_ags = {ag.name for ag in agents}
        if len(uniq_ags) < 2:
            continue
        if uniq_ags < set(channels):
            continue
        new_stmts.append(stmt)
    return new_stmts


def print_statistics(statements):
    counts = sorted([(k, len(s[0])) for k, s in statements.items()],
                    key=lambda x: x[1], reverse=True)
    raw_counts = [c[1] for c in counts]
    missing = [c[0] for c in counts if c[1] == 0]
    print(f'No statements for channels: {", ".join(missing)}')
    print(f'{counts[0][1]} statements for the top channel {counts[0][0]}')
    print(f'{numpy.mean(raw_counts)} statements on average per channel')


if __name__ == '__main__':
    fname = 'data/gene_list.txt'
    with open(fname, 'r') as fh:
        channels = [l.strip() for l in fh.readlines()]

    fname = 'data/IDG_target_final.csv'
    df = pandas.read_csv(fname)
    print('Read a total of %d rows from %s' % (len(df), fname))
    df = df[df['idgFamily'] == 'Ion Channel']
    print('Filtered to %d ion channels' % len(df))
    all_channel_gene_names = sorted(list(df['gene']))
    print('Read a total of %d channels from %s' % (len(channels), fname))

    stmts_with_counts = get_channel_statements(channels)
    with open('ion_channel_stmts_v1.pkl', 'wb') as fh:
        pickle.dump(stmts_with_counts, fh)
    print_statistics(stmts_with_counts)
