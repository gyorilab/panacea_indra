import os
import pickle
from collections import defaultdict
from indra.sources import tas
from indra.statements import *
from indra.databases import get_identifiers_url
from indra.assemblers.html import HtmlAssembler


with open('../ion_channel_regulators.pkl', 'rb') as fh:
    neg_regs, non_neg_regs = pickle.load(fh)


with open('../ion_channel_stmts_v4.pkl', 'rb') as fh:
    stmts_by_channel = pickle.load(fh)


def get_key(agent):
    gr = agent.get_grounding()
    if gr == (None, None):
        gr = ('NAME', agent.name)
    return gr


def find_reg_non_regs(reg_channels, non_reg_channels):
    # Take the intersection of known regulators of channels we want to target
    joint_regs = set(neg_regs[reg_channels[0]])
    for channel in reg_channels:
        joint_regs &= set(neg_regs[channel].keys())
    # Remove the union of known regulators of channels we don't want to target
    remove_regs = set()
    for channel in non_reg_channels:
        remove_regs |= set(neg_regs[channel].keys())
    keys = joint_regs - remove_regs
    # Recover the Agents behind the keys
    agents = [neg_regs[reg_channels[0]][k] for k in keys]
    return agents


def get_agent_urls(agent):
    urls = [(db, get_identifiers_url(db, id))
            for db, id in agent.db_refs.items() if db != 'TEXT']
    urls = [u for u in urls if u[1] is not None]
    return urls


def get_reg_stmts(reg_agent, channel):
    reg_stmts = []
    for stmt in stmts_by_channel[channel]:
        if isinstance(stmt, (Inhibition, DecreaseAmount)):
            if stmt.subj.name == reg_agent.name:
                reg_stmts.append(stmt)
    return reg_stmts


def assemble_html(stmts, fname_key):
    ha = HtmlAssembler(stmts)
    ha.make_model()
    ha.save_model('%s.html' % fname_key)


if __name__ == '__main__':
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
