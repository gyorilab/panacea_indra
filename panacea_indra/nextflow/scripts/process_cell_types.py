import os
import sys
import pickle
import logging
import pandas as pd
import enzyme_client
from collections import defaultdict
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name


logger = logging.getLogger('receptor_ligand_interactions')

# Functions

def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def process_seurat_csv(infile, fc):
    """ Process Seurat dataframe and only filter in
    genes with the given Fold change """
    l_df = pd.read_csv(infile, header=0, sep=",")
    l_df.columns = l_df.columns.str.replace('Unnamed: 0', 'Genes')
    filtered_df = l_df[l_df['avg_logFC'] > 0.25][['Genes', 'avg_logFC']]
    filtered_df = filtered_df.sort_values(by='avg_logFC', ascending=False)
    filtered_dict = {}
    for r, c in filtered_df.iterrows():
        filtered_dict[c[1]] = c[0]
    # Volcano plot of DE genes
    #_plot_de_genes(l_df, output_dir)
    # return set(filtered_markers)
    return filtered_dict


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


def make_rnk(infile):
    df = pd.read_csv(infile, header=0, sep=",")
    df.columns = df.columns.str.replace('Unnamed: 0', 'Genes')
    df = df.loc[0:, ['Genes', 'p_val']]
    return df


def get_de_product_list(de_enzyme_product_list,
                        de_enzyme_stmts):
    if len(de_enzyme_product_list) > 1:
        de_enzyme_product_list = pd.merge(de_enzyme_stmts, de_enzyme_product_list,
                                          on=['Enzyme', 'Interaction', 'product', 'logFC'],
                                          how="outer").fillna('')
        return de_enzyme_product_list.sort_values(by='logFC', ascending=False)

    elif len(de_enzyme_product_list) < 1:
        de_enzyme_product_list = de_enzyme_stmts
        return de_enzyme_product_list


def get_pain_mol():
    PAIN_SIGNAL_MOL = {
        "Prostaglandins": "CHEBI:26333",
        "Brandykinin": "CHEBI:3165"
    }

    CHEBI_LIST = {}
    CHEBI_NAMES = {}
    for compounds, chebi_id in PAIN_SIGNAL_MOL.items():
        CHEBI_LIST[compounds] = \
            [children[1] for children in
             bio_ontology.get_children('CHEBI',
                                       chebi_id)]

        CHEBI_NAMES[compounds] = \
            [bio_ontology.get_name('CHEBI', ids)
             for ids in CHEBI_LIST[compounds]]

    return CHEBI_NAMES


def get_hashes_by_gene_pair(df, ligand_genes, receptor_genes):
    hashes_by_gene_pair = defaultdict(set)
    l_genes = list(ligand_genes.values())
    l_logFC = list(ligand_genes.keys())

    for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
        if a in l_genes and b in receptor_genes:
            logFC = l_logFC[l_genes.index(a)]
            hashes_by_gene_pair[(a, b, logFC)].add(hs)
    return hashes_by_gene_pair





if __name__ == '__main__':
    INPUT = sys.argv[1]
    OUTPUT = sys.argv[2]
    ALL_ENZYMES = sys.argv[3]
    FULL_LIGAND_SET = sys.argv[4]
    TARGETS_BY_DRUG = sys.argv[5]

    PAIN_MOL_NAMES = get_pain_mol()

    # Read channels from nextflow
    with open(ALL_ENZYMES, 'rb') as fh:
        all_enzymes = pickle.load(fh)

    with open(FULL_LIGAND_SET, 'rb') as fh:
        full_ligand_set = pickle.load(fh)

    with open(os.path.join(OUTPUT, 'receptors_in_data.pkl'), 'rb') as fh:
        receptors_in_data = pickle.load(fh)

    with open(TARGETS_BY_DRUG, 'rb') as fh:
        targets_by_drug = pickle.load(fh)



    mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}

    #IMMUNE_CELLTYPE_LIST = ['DCs',
    #                        'Dermal Macs']



    IMMUNE_CELLTYPE_LIST = ['DCs',
                            'Dermal Macs',
                            'M2a',
                            'M2b',
                            'Monocytes',
                            'Resident Mac',
                            'Mast cells'
                            ]


    ligands_df = pd.DataFrame(columns=['Genes', 'p_val'])
    ligands_FC = {}
    enzymes_FC = {}
    cell_type_markers = {}
    possible_drug_targets = set()
    de_enzyme_list = set()
    de_enzyme_product_list = set()
    possible_en_drug_targets = defaultdict(set)
    

    for cell_type in IMMUNE_CELLTYPE_LIST:
        logger.info('Processing %s' % cell_type)
        enzyme_product_dict = defaultdict(set)

        INDRA_DB_PKL = os.path.join(INPUT, 'db_dump_df.pkl')
        # Load the INDRA DB DF
        df = load_indra_df(INDRA_DB_PKL)

        output_dir = os.path.join(OUTPUT, cell_type)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # read the input (immune cell type) ligand file
        cell_type_full = 'TwoGroups_DEG1_%s_AJ' % cell_type
        LIGANDS_INFILE = os.path.join(INPUT, '%s.csv' % cell_type_full)
        # Extract markers from seurat dataframe with logFC >= 0.25 and
            # pval <= 0.05
        seurat_ligand_genes = process_seurat_csv(LIGANDS_INFILE,
                                                fc=0.25)

        # Pool all ligands along with its respective logFC and create a rank file
        # for GSEA
        ligands_df = pd.concat([ligands_df, make_rnk(LIGANDS_INFILE)],
                               ignore_index=True)

        if len(seurat_ligand_genes) == 0:
                logger.info('Skipping %s' % cell_type)
                continue

        # Get logFC as key and ligands as values
        ligand_genes = ligand_mgi_to_hgnc_name(seurat_ligand_genes)

        # Retain only ligands
        ligands_in_data = {k: l for k, v in ligand_genes.items() for l in v
                           if l in full_ligand_set}

        # Retain only enzymes
        enzymes_in_data = {k: next(iter(v)) for k, v in ligand_genes.items()
                           if next(iter(v)) in all_enzymes}

        # Keep all ligands with FC from all cell types
        ligands_FC.update(ligands_in_data)
        # Keep all enzymes with FC from all cell typesa
        enzymes_FC.update(enzymes_in_data)

        cell_type_markers[cell_type] = ligands_in_data
        cell_type_markers[cell_type].update(enzymes_in_data)

        # Get enzyme products by taking pathway commons DB as reference
        de_enzyme_product_interaction = enzyme_client.get_enzyme_products(enzymes_in_data)


        with open(cell_type+'_de_enzyme_product_interaction.pkl', 'wb') as fh:
            pickle.dump(de_enzyme_product_interaction, fh)

        # write de enzyme products to a file
        de_enzyme_product_interaction.to_csv(os.path.join(output_dir, cell_type + "_de_enzymes_product_interaction.tsv"),
                               sep="\t", header=True, index=False)

        # Store enzymes as keys and its respective product
        # as values
        for r, c in de_enzyme_product_interaction.iterrows():
            enzyme_product_dict[(c[0])].add(c[2])

        # Keep merging enzyme products interactions from all the celltypes
        de_enzyme_product_list = get_de_product_list(de_enzyme_product_list,
                                                     de_enzyme_product_interaction)

        # Get enzyme interactions with pain molecules for the cell type
        pain_interactions = enzyme_client.get_pain_interactions(de_enzyme_product_interaction,
                                                                PAIN_MOL_NAMES)

        with open(cell_type+'_pain_interactions.pkl', 'wb') as fh:
            pickle.dump(pain_interactions, fh)

        # write de enzyme pain interactions to a file
        pain_interactions.to_csv(os.path.join(output_dir, cell_type + "_enzyme_pain_interactions.tsv"),
                                 sep="\t", header=True, index=False)


        # Get the union of all the enzymes in the data
        # and include them in drug target set
        possible_drug_targets |= set(enzymes_in_data.values())


        # Get the union of all the enzymes in the data
        # from all the cell types
        de_enzyme_list |= set(enzymes_in_data.values())

        # Make a data frame of ligands and logFC
        # for the respective cell type
        ligands_fc_df = {'Ligands': [*ligands_in_data.values()],
                         'logFC': [*ligands_in_data.keys()]}
        ligands_fc_df = pd.DataFrame(ligands_fc_df)

                # Save the dataframe into a csv
        ligands_fc_df.to_csv(os.path.join(output_dir, cell_type + "_ligands_fc.csv"),
                             header=True, index=False)
        logger.info(f'Loaded {len(ligands_in_data)} ligand genes from data')
        logger.info(f'Loaded {len(receptors_in_data)} receptor genes from data')

        # Now get INDRA DB Statements for the receptor-ligand pairs
        hashes_by_gene_pair = get_hashes_by_gene_pair(df, ligands_in_data,
                                                      receptors_in_data)

        with open(cell_type+'_hashes_by_gene_pair.pkl', 'wb') as fh:
            pickle.dump(hashes_by_gene_pair, fh)

        with open('de_enzyme_product_list.pkl', 'wb') as fh:
            pickle.dump(de_enzyme_product_list, fh)

        # Create a dictionary of enzymes and its foldchange
        for fc, en in enzymes_in_data.items():
            possible_en_drug_targets[(fc)].add(en)



        # Get the list of ligands from the cell type
        #ligands_in_data = list(ligands_fc_df['Ligands'])

        with open(cell_type+'_ligands_in_data.pkl', 'wb') as fh:
            pickle.dump(ligands_in_data, fh)

    
    with open(os.path.join(OUTPUT, 'enzymes_in_data.pkl'), 'wb') as fh:
        pickle.dump(enzymes_in_data, fh)

    with open('enzyme_product_dict.pkl', 'wb') as fh:
        pickle.dump(enzyme_product_dict, fh)

    with open('possible_en_drug_targets.pkl', 'wb') as fh:
        pickle.dump(possible_en_drug_targets, fh)

    with open('enzyme_possible_drug_targets.pkl', 'wb') as fh:
        pickle.dump(possible_drug_targets, fh)

    with open('ligands_FC.pkl', 'wb') as fh:
        pickle.dump(ligands_FC, fh)

    with open('enzymes_FC.pkl', 'wb') as fh:
        pickle.dump(enzymes_FC, fh)














