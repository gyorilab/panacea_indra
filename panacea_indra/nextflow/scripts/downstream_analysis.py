import os
import re
import sys
import csv
import json
import pickle
import datetime
import argparse
import openpyxl
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from indra.databases.uniprot_client import um
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name


mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}


def fix_dates(gene_names):
    replacements = {
        datetime.datetime(2020, 3, 7, 0, 0): 'March7',
        datetime.datetime(2020, 3, 2, 0, 0): 'March2',
        datetime.datetime(2020, 3, 4, 0, 0): 'March4',
        datetime.datetime(2020, 3, 5, 0, 0): 'March5',
        datetime.datetime(2020, 3, 6, 0, 0): 'March6',
        datetime.datetime(2020, 3, 9, 0, 0): 'March9',
        datetime.datetime(2020, 3, 8, 0, 0): 'March8',
        datetime.datetime(2020, 3, 11, 0, 0): 'Mar11',
        datetime.datetime(2020, 9, 1, 0, 0): 'Sept1',
        datetime.datetime(2020, 9, 2, 0, 0): 'Sept2',
        datetime.datetime(2020, 9, 3, 0, 0): 'Sept3',
        datetime.datetime(2020, 9, 4, 0, 0): 'Sept4',
        datetime.datetime(2020, 9, 5, 0, 0): 'Sept5',
        datetime.datetime(2020, 9, 6, 0, 0): 'Sept6',
        datetime.datetime(2020, 9, 7, 0, 0): 'Sept7',
        datetime.datetime(2020, 9, 8, 0, 0): 'Sept8',
        datetime.datetime(2020, 9, 9, 0, 0): 'Sept9',
        datetime.datetime(2020, 9, 10, 0, 0): 'Sept10',
        datetime.datetime(2020, 9, 11, 0, 0): 'Sept11',
        datetime.datetime(2020, 9, 15, 0, 0): 'Sept15',
    }
    fixed_gene_names = set()
    for gene_name in gene_names:
        if isinstance(gene_name, datetime.datetime):
            fixed_gene_names.add(replacements[gene_name])
        else:
            fixed_gene_names.add(gene_name)
    return fixed_gene_names


def get_pain_phenotype(lg, pain_db):
    r_phenotype = defaultdict(set)
    for r, c in pain_db.iterrows():
        if isinstance(pain_db.iloc[r]['gene_symbols'], str):
            l = set(pain_db.iloc[r]['gene_symbols'].split(","))
            pheno = pain_db.iloc[r]['phenotype_description']
            for g in l:
                if g in lg:
                    r_phenotype[(g)].add(pheno)
    return r_phenotype


def read_workbook(workbook):
    """ This function takes Excel workbook as an input and
    returns ligand and receptor gene list respectively.
    Input: Excel workbook with single(2 columns) or two sheets
    Condition: considers first column/sheet as ligand genes and second
    column/shet as receptor genes
    """
    ligands_sheet = 'updated list of ligands '
    receptors_sheet = 'RPKM > 1.5 cfiber'
    wb = openpyxl.load_workbook(workbook)
    ligands = fix_dates(set([row[0].value for row in wb[ligands_sheet]][1:]))
    receptors = fix_dates(set([row[0].value
                               for row in wb[receptors_sheet]][1:]))
    return ligands, receptors


def ligand_mgi_to_hgnc_name(seurat_ligand_genes):
    filtered_mgi = defaultdict(set)
    for logfc, gene in seurat_ligand_genes.items():
        if gene in mouse_gene_name_to_mgi:
            filtered_mgi[(gene, logfc)].add(mouse_gene_name_to_mgi[gene])

    hgnc_gene_dict = defaultdict(set)
    seen_genes = set()
    for key, value in filtered_mgi.items():
        mgi_id = next(iter(value))
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_symbol = get_hgnc_name(hgnc_id)
        if hgnc_symbol not in seen_genes:
            hgnc_gene_dict[(key[1])].add(hgnc_symbol)
        else:
            pass
        seen_genes.add(hgnc_symbol)
    return hgnc_gene_dict


def mgi_to_hgnc_name(gene_list):
    """Convert given mouse gene symbols to HGNC equivalent symbols"""
    filtered_mgi = {mouse_gene_name_to_mgi[gene] for gene in gene_list
                    if gene in mouse_gene_name_to_mgi}
    hgnc_gene_set = set()
    for mgi_id in filtered_mgi:
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_gene_set.add(get_hgnc_name(hgnc_id))
    return hgnc_gene_set


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--human_pain_db")
    parser.add_argument("--receptor_genes_go")
    parser.add_argument("--targets_by_drug")
    args = parser.parse_args()
    
    INPUT = args.input
    OUTPUT = args.output
    HUMAN_PAIN_DB = args.human_pain_db
    RECEPTOR_GENES_GO = args.receptor_genes_go
    TARGETS_BY_DRUG = args.targets_by_drug
    DATA_SPREADSHEET = os.path.join(INPUT, 'Neuroimmune gene list .xlsx')

    with open(os.path.join(OUTPUT, 'receptors_in_data.pkl'), 'rb') as fh:
        receptors_in_data = pickle.load(fh)

    with open(RECEPTOR_GENES_GO, 'rb') as fh:
        receptor_genes_go = pickle.load(fh)

    with open(os.path.join(INPUT, 'db_dump_df.pkl'), 'rb') as fh:
        indra_sif = pickle.load(fh)

    with open(TARGETS_BY_DRUG, 'rb') as fh:
        targets_by_drug = pickle.load(fh)


    _, raw_receptor_genes = read_workbook(DATA_SPREADSHEET)
    receptor_genes = mgi_to_hgnc_name(raw_receptor_genes)

    # Read human pain DB
    human_pain_df = pd.read_csv(HUMAN_PAIN_DB, sep="\t", header=0)

    # Create a dictionary of receptors as keys
    # and its respective pain phenotypes as values
    r_phenotype_dict = get_pain_phenotype(receptors_in_data,
	                                      human_pain_df)


    ### Downstream analysis

    # Human pain genetics db
    gene_symbol = [str(g).split(",") for g in human_pain_df['gene_symbols']]
    pain_phenotype = [str(p).split(",") for p in human_pain_df['phenotype_description']]
    pheno_direction = human_pain_df['phenotype_direction']
    pain_pheno_dict = defaultdict(set)

    # Import ctd database
    # Opening JSON file and loading the data
    # into the variable data
    with open(os.path.join(INPUT, 'ctd_db.json')) as json_file:
        data = json.load(json_file)

    ctd_db_json = data['associations']


    # now we will open a file for writing
    data_file = open(os.path.join(INPUT, 'ctd_db.csv'), 'w')

    # create the csv writer object
    csv_writer = csv.writer(data_file)

    # Counter variable used for writing
    # headers to the CSV file
    count = 0

    for emp in ctd_db_json:
        if count == 0:
            # Writing headers of CSV file
            header = emp.keys()
            csv_writer.writerow(header)
            count += 1

        # Writing data of CSV file
        csv_writer.writerow(emp.values())
    data_file.close()


    # Read the csv back and clean it
    ctd_csv = pd.read_csv(os.path.join(INPUT, 'ctd_db.csv'), header=0)
    count = 0
    # Clean the gene column
    for r, c in ctd_csv.iterrows():
        gene = re.search("symbol': '([A-Z0-9]+)", c[0]).group(1)
        ctd_csv['gene'][count] = gene
        count += 1


    # Creating a dictionary of ctd genes
    # as keys and its respective Standarized values
    # as values
    ctd_geneset = dict()
    for r, c in ctd_csv.iterrows():
        ctd_geneset[c[0]] = c[2]

    downstream_hits = defaultdict(set)
    downstream_evidence = dict()
    downstream_stmts = defaultdict(set)
    rc_in_data_dwnstrm = set()
    upstream = defaultdict(set)
    not_upstream = defaultdict(set)


    for a, b, ev, hs in zip(indra_sif.agA_name,
                            indra_sif.agB_name,
                            indra_sif.evidence_count,
                            indra_sif.stmt_hash):
        upstream[(b)].add(a)
        if a in receptor_genes_go:
            downstream_hits[(b)].add(a)
            if b not in downstream_evidence:
                downstream_evidence[b] = ev
            else:
                downstream_evidence[b] = downstream_evidence[b] + ev
            downstream_stmts[(b)].add(hs)

        if a in receptors_in_data:
            rc_in_data_dwnstrm.add(b)


    targets_keys = list(targets_by_drug.keys())
    targets_values = list(targets_by_drug.values())

    ds_in_data = defaultdict(set)
    drug_names = defaultdict(set)
    drug_no_names = defaultdict(set)

    for drug, targets in targets_by_drug.items():
        for t in targets:
            if t in downstream_hits.keys():
                if drug[0].startswith('CHEMBL'):
                    drug_no_names[(t)].add(drug[0])
                else:
                    drug_names[(t)].add(drug[0])


    # Below code is to aggregate all the
    # downstream hits and its pain pheno direction
    ds_direction_dict = dict()
    pheno_direction = human_pain_df['phenotype_direction']
    for k in downstream_hits.keys():
        for g, p, pdir in zip(gene_symbol,
                              pain_phenotype,
                              pheno_direction):
            if k in g:
                p = "".join(p)
                if k + ":" + p + ":" + pdir in ds_direction_dict:
                    ds_direction_dict[k + ":" + p + ":" + pdir] = ds_direction_dict[k + ":" + p + ":" + pdir] + 1
                else:
                    ds_direction_dict[k + ":" + p + ":" + pdir] = 1
    d = defaultdict(set)

    for t in ds_direction_dict:
        d[re.match("([a-zA-Z0-9-]+)(.*)", t).group(1)].add(
            re.match("([a-zA-Z0-9-]+:)(.*)", t).group(2) + ":" + str(ds_direction_dict[t]))

    # Creating a dictionary of downstream targets as keys
    # its values as pain phenotype description
    for k in downstream_hits.keys():
        for g, p in zip(gene_symbol, pain_phenotype):
            if k in g:
                pain_pheno_dict[(k)].add("".join(p))


    downstream_df = []
    for ds, us in downstream_hits.items():
        if ds in receptor_genes:
            total_interactome = len(receptors_in_data)
            total_non_interactome = len(receptor_genes_go - receptors_in_data)
            interactome_genes = len(receptors_in_data & downstream_hits[ds])
            non_interactome_genes = len(downstream_hits[ds] - receptors_in_data)

            pval = fisher_exact([[interactome_genes, non_interactome_genes],
                                 [total_interactome, total_non_interactome]],
                                alternative='greater')[1]

            if ds in drug_names:
                drugs_names = ", ".join(drug_names[ds])
            else:
                drugs_names = 'NA'

            if ds in drug_no_names:
                drugs_no_names = ", ".join(drug_no_names[ds])
            else:
                drugs_no_names = 'NA'

            if ds in ctd_geneset:
                ctd_gene = 1
                std_val = ctd_geneset[ds]
            else:
                ctd_gene = 0
                std_val = 0

            if ds in d:
                pain_phenotype_des = ", ".join(d[ds])
            else:
                pain_phenotype_des = 'NA'

            downstream_df.append(
                {
                    "Upstream genes": ", ".join(us),
                    "Upstream in interacome": ", ".join(receptors_in_data & downstream_hits[ds]),
                    "Downstream target": ds,
                    "CTD_DB": ctd_gene,
                    "CTD_Std_Val": std_val,
                    "Pain description": pain_phenotype_des,
                    "Named drugs": drugs_names,
                    "Unnamed drugs": drugs_no_names,
                    "evidence count": downstream_evidence[ds],
                    "Total upstream": len(us),
                    "Total upstream in interactome": interactome_genes,
                    "pval": pval,
                    "Statement hash": ', '.join(str(x) for x in downstream_stmts[ds])
                }
            )
    downstream_df = pd.DataFrame(downstream_df).sort_values(by=['pval'],
                                                            ascending=True)

    downstream_df.to_csv(os.path.join(OUTPUT, "downstream_hits_filtered_pval.csv"),
                         sep=",", header=True, index=False)
