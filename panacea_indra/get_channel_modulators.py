import os
import pickle
from collections import defaultdict
from indra.sources import tas
from indra.statements import *
from indra.databases import get_identifiers_url


with open('../ion_channel_regulators.pkl', 'rb') as fh:
    neg_regs, non_neg_regs = pickle.load(fh)


def get_key(agent):
    gr = agent.get_grounding()
    if gr == (None, None):
        gr = ('NAME', agent.name)
    return gr


def find_reg_non_regs(reg_channels, non_reg_channels):
    joint_regs = set(neg_regs[reg_channels[0]])
    for channel in reg_channels:
        joint_regs &= set(neg_regs[channel].keys())
    remove_regs = set()
    for channel in non_reg_channels:
        remove_regs |= set(neg_regs[channel].keys())
    keys = joint_regs - remove_regs
    agents = [neg_regs[reg_channels[0]][k] for k in keys]
    return agents


def get_agent_urls(agent):
    urls = [get_identifiers_url(db, id) for db, id in agent.db_refs.items()
            if db != 'TEXT']
    urls = [u for u in urls if u is not None]
    return urls


if __name__ == '__main__':
    with open('../ion_channel_stmts_v3.pkl', 'rb') as fh:
        stmts_by_channel = pickle.load(fh)
    tp = tas.process_from_web(affinity_class_limit=10)

    neg_regs = defaultdict(dict)
    non_neg_regs = defaultdict(dict)
    for channel, (stmts, _, _) in stmts_by_channel.items():
        stmts = [s for s in stmts if
                 isinstance(s, (Inhibition, DecreaseAmount))]
        stmts = [s for s in stmts if s.obj.name == channel]
        for stmt in stmts:
            neg_regs[channel][get_key(stmt.subj)] = stmt.subj

    for stmt in tp.statements:
        if stmt.obj.name in neg_regs:
            if stmt.evidence[0].annotations['class_min'] in \
                    {'100nM < Kd < 1uM', 'Kd < 100nM'}:
                neg_regs[stmt.obj.name][get_key(stmt.subj)] = stmt.subj
            else:
                non_neg_regs[stmt.obj.name][get_key(stmt.subj)] = stmt.subj

    with open('../ion_channel_regulators.pkl', 'wb') as fh:
        pickle.dump((neg_regs, non_neg_regs), fh)
