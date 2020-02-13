import gilda
import pickle
from collections import defaultdict
from indra.sources import tas
from indra.statements import *


def get_subj_agents(stmts):
    others = set()
    for stmt in stmts:
        others.add(get_key(stmt.subj))
    return others


def get_key(agent):
    gr = agent.get_grounding()
    if gr == (None, None):
        gr = ('NAME', agent.name)
    return gr


def find_joint_regs(channels):
    joint_regs = set(neg_regs[channels[0]])
    for channel in channels:
        joint_regs &= set(neg_regs[channel])
    return joint_regs


def find_reg_non_regs(reg_channels, non_reg_channels):
    joint_regs = find_joint_regs(reg_channels)
    remove_regs = set()
    for channel in non_reg_channels:
        remove_regs |= set(neg_regs[channel])
    return joint_regs - remove_regs


if __name__ == '__main__':
    with open('ion_channel_stmts_v3.pkl', 'rb') as fh:
        stmts_by_channel = pickle.load(fh)
    tp = tas.process_from_web(affinity_class_limit=10)

    neg_regs = {}
    non_neg_regs = {}
    for channel, (stmts, _, _) in stmts_by_channel.items():
        stmts = [s for s in stmts if
                 isinstance(s, (Inhibition, DecreaseAmount))]
        stmts = [s for s in stmts if s.obj.name == channel]
        neg_regs[channel] = list(get_subj_agents(stmts))
        non_neg_regs[channel] = []

    for stmt in tp.statements:
        if stmt.obj.name in neg_regs:
            if stmt.evidence[0].annotations['class_min'] in \
                    {'100nM < Kd < 1uM', 'Kd < 100nM'}:
                neg_regs[stmt.obj.name].append(get_key(stmt.subj))
            else:
                non_neg_regs[stmt.obj.name].append(get_key(stmt.subj))

    with open('ion_channel_regulators.pkl', 'wb') as fh:
        pickle.dump((neg_regs, non_neg_regs), fh)