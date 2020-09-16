import re
import requests
import pandas as pd
from indra.ontology.bio import bio_ontology


ENZYME_URL = 'ftp://ftp.expasy.org/databases/enzyme/enzyme.dat'
PC_SIF_URL = ('https://www.pathwaycommons.org/archives/PC2/v12/'
              'PathwayCommons12.Detailed.hgnc.sif.gz')


def load_enzyme_data(dat_str):
    enzyme_data = {}
    enzyme_class = None
    enzyme_entries = []
    for line in dat_str.split('\n'):
        g = re.match(r'^ID\s+([0-9\.]+)$', line)
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
    human_data = {}
    for enz_class, enz_entries in enzyme_data.items():
        human_entries = [e for e in enz_entries if e[1].endswith('_HUMAN')]
        human_data[enz_class] = human_entries
    return human_data


if __name__ == "__main__":
    filtered_stmts = []
    df = pd.read_csv(PC_SIF_URL, sep='\t', header=None)
    df = df[df[1] == 'controls-production-of']

    pain_signal_mol = {
        "Prostaglandins": "CHEBI:26333",
        "Brandykinin": "CHEBI:3165"
    }

    chebi_list = {}
    for compounds, chebi_id in pain_signal_mol.items():
        chebi_list[compounds] = [children[1] for children in
                                 bio_ontology.get_children('CHEBI',
                                                           chebi_id)]
    df = df[df[2].isin(chebi_list)]
    chebi_stmts = [
        {'Enzyme': row[0],
         'Statement': row[1],
         'CHEBI_ID': row[2],
         'CHEBI_Name': bio_ontology.get_name('CHEBI', row[2])}
        for _, row in df.iterrows()]
    df = pd.DataFrame(chebi_stmts)
    df.to_csv("enzyme_interactions.tsv", sep="\t",
              header=True, index=False)
