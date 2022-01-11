import os
import sys
import pickle
import pandas as pd


def get_small_mol_report(targets_by_drug, potential_targets, fname):
    df = []
    for drug, targets in targets_by_drug.items():
        targets_in_data = targets & potential_targets
        if not targets_in_data:
            continue
        df.append(
            {
                "Drug": drug[0],
                "ID": '%s:%s' % (drug[1]),
                "Named": 0 if drug[0].startswith('CHEMBL') else 1,
                "Score": "{:.3f}".format(len(targets_in_data) / len(targets)),
                "Number of targets in data": len(targets_in_data),
                "Targets in data": ", ".join(sorted(targets_in_data)),
                "Other targets": ", ".join(sorted(targets - targets_in_data)),
            }
        )
    df = pd.DataFrame(df).sort_values(by=['Score', 'Number of targets in data',
                                          'Named'],
                                      ascending=False)
    df.to_csv(fname, sep="\t", header=True, index=False)
    return df


if __name__ == '__main__':

	INPUT = sys.argv[1]
	OUTPUT = sys.argv[2]
	TARGETS_BY_DRUG = sys.argv[3]
	POSSIBLE_DRUG_TARGETS = sys.argv[4]
	POSSIBLE_DB_DRUG_TARGETS = sys.argv[5]

	with open(TARGETS_BY_DRUG, 'rb') as fh:
		targets_by_drug = pickle.load(fh)

	with open(POSSIBLE_DRUG_TARGETS, 'rb') as fh:
		possible_drug_targets = pickle.load(fh)

	with open(POSSIBLE_DB_DRUG_TARGETS, 'rb') as fh:
		possible_db_drug_targets = pickle.load(fh)


	get_small_mol_report(targets_by_drug, possible_drug_targets,
		os.path.join(OUTPUT, 'drug_targets.tsv'))

	get_small_mol_report(targets_by_drug, possible_db_drug_targets,
		os.path.join(OUTPUT, 'drug_targets_db.tsv'))