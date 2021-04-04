import os
import re
import sys
import pickle
import argparse
import pandas as pd
from collections import defaultdict, OrderedDict





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--ligands_by_receptor", nargs='+')
    parser.add_argument("--possible_en_drug_targets")
    parser.add_argument("--targets_by_drug")
    parser.add_argument("--enzyme_product_dict")
    args = parser.parse_args()


    INPUT = args.input
    OUTPUT = args.output
    LIGANDS_BY_RECEPTOR = args.ligands_by_receptor
    POSSIBLE_EN_DRUG_TARGETS = args.possible_en_drug_targets
    TARGETS_BY_DRUG = args.targets_by_drug
    ENZYME_PRODUCT_DICT = args.enzyme_product_dict




    ### Enzyme product drug interaction
    ligands_by_receptor = {}

    for lg_by_rg in LIGANDS_BY_RECEPTOR:
        try:
            with open(lg_by_rg, 'rb') as fh:
                rg_lg = pickle.load(fh)

            for k, v in rg_lg.items():
                if k not in ligands_by_receptor:
                    ligands_by_receptor[k] = v
                else:
                    ligands_by_receptor[k].update(v)
        except EOFError:
            continue
    



    with open(POSSIBLE_EN_DRUG_TARGETS, 'rb') as fh:
        possible_en_drug_targets = pickle.load(fh)


    with open(TARGETS_BY_DRUG, 'rb') as fh:
        targets_by_drug = pickle.load(fh)

    with open('enzyme_product_dict.pkl', 'rb') as fh:
        enzyme_product_dict = pickle.load(fh)

            
    # full drug interactions
    possible_en_drug_targets = OrderedDict(sorted(possible_en_drug_targets.items(),
                                                  reverse=True))
    # Load logFC of enzyme drug targets into a list
    fc_en_drug_targets = list(possible_en_drug_targets.keys())

    # Load enzymes into list
    en_drug_targets = list(possible_en_drug_targets.values())

    en_targets_in_data = defaultdict(set)
    rc_targets_in_data = defaultdict(set)
    fc_in_data = defaultdict(set)

            
    for drug, targets in targets_by_drug.items():
        for t in targets:
            for e in en_drug_targets:
                if t in e:
                    fc = fc_en_drug_targets[en_drug_targets.index(e)]
                    fc_in_data[(drug[0])].add(fc)
                    en_targets_in_data[(drug[0])].add(t)
                if t in ligands_by_receptor.keys():
                    rc_targets_in_data[(drug[0])].add(t)




    drug_interaction_list = []
    for drug, targets in targets_by_drug.items():
        intermediates = set()
        l = []
        if drug[0] in en_targets_in_data.keys():
            en_target = list(en_targets_in_data[drug[0]])
            for enzymes in en_target:
                intermediates.update(enzyme_product_dict[enzymes])
            intermediates = list(intermediates)

            for val in intermediates:
                if val != None:
                    l.append(val)
            avgFC = sum(fc_in_data[drug[0]]) / len(fc_in_data[drug[0]])
        else:
            en_target = ''
            avgFC = 0

        if drug[0] in rc_targets_in_data.keys():
            rc_target = rc_targets_in_data[drug[0]]
        else:
            rc_target = ""

        other_targets = targets - set(en_target) - set(rc_target)

        drug_interaction_list.append(
            {
                "Drug": drug[0],
                "Named": 0 if drug[0].startswith('CHEMBL') else 1,
                "Enzyme targets": ", ".join(en_target),
                "Enzyme products": ", ".join(l),
                "Receptor targets": ", ".join(rc_target),
                "Score": "{:.3f}".format((len(en_target) + len(rc_target)) / len(targets)),
                "Other targets": ", ".join(sorted(other_targets)),
                "Other target hits": "{:.3f}".format(len(other_targets)),
                "Total hits": "{:.3f}".format((len(en_target) + len(rc_target) + len(other_targets))),
                "avgFC": float(avgFC)
            }
        )
    pd.DataFrame(drug_interaction_list).sort_values(by=['avgFC'], ascending=False).to_csv(os.path.join(
        OUTPUT, "ranked_enzyme_drug_target.csv"))
