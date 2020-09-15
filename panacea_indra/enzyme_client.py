import os
import re
import logging
import collections
import pandas as pd
from indra.ontology.bio import bio_ontology


ENZYME_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           './data/enzyme.dat')
PATHWAY_HGNC_SIF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           './data/PathwayCommons12.Detailed.hgnc.sif')


def load_enzyme_data(filename):
    with open(enzyme_file) as f:
        lines = f.readlines()
    enzyme_data = collections.OrderedDict()
    enzyme_class = None
    enzyme_entries = []
    for line in lines:
        g = re.match('^ID\s+([0-9\.]+)$', line)
        if g:
            enzyme_class = g.groups()[0]
            enzyme_entries = []
            continue
        if line[0:2] == '//':
            enzyme_data[enzyme_class] = enzyme_entries
            enzyme_class = None
            continue
        if enzyme_class is not None and line[0:2] == 'DR':
            entries = line[2:].strip().split(';')
            for entry in entries:
                if entry:
                    (up_id, up_mnemonic) = entry.strip().split(',')
                    up_id = up_id.strip()
                    up_mnemonic = up_mnemonic.strip()
                    enzyme_entries.append((up_id, up_mnemonic))
    return enzyme_data


def filter_human(enzyme_data):
    human_data = collections.OrderedDict()
    for enz_class, enz_entries in enzyme_data.items():
        human_entries = [e for e in enz_entries if e[1].endswith('_HUMAN')]
        human_data[enz_class] = human_entries
    return human_data


if __name__ == "__main__":
    filtered_stmts = []
    with open(PATHWAY_HGNC_SIF) as fh:
        for lines in fh:
            if re.search("controls-production-of", lines):
                filtered_stmts.append(lines.rstrip())
    len(filtered_stmts)

    pain_signal_mol = {
        "Prostaglandins": "CHEBI:26333",
        "Brandykinin": "CHEBI:3165"
    }

    chebi_list = collections.OrderedDict()
    for compounds, chebi_id in pain_signal_mol.items():
        chebi_list[compounds] = [children[1] for children in
                                 bio_ontology.get_children('CHEBI',
                                                           chebi_id)]

    chebi_stmts = []
    for index, chebi_entries in chebi_list.items():
        for entry in chebi_entries:
            for stmts in filtered_stmts:
                if re.search(entry, stmts):
                    splitted_values = stmts.split("\t")
                    chebi_name = bio_ontology.get_name('CHEBI',
                                                       splitted_values[2])
                    chebi_stmts.append(
                        {
                            'Enzymes': splitted_values[0],
                            'Statement': splitted_values[1],
                            'CHEBI_ID': splitted_values[2],
                            'CHEBI_Name': chebi_name
                        }
                    )
    df = pd.DataFrame(chebi_stmts)
    df.to_csv("enzyme_interactions.tsv", sep="\t",
              header=True, index=False)