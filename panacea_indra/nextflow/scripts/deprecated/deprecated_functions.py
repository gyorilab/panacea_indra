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

PAIN_MOL_NAMES = get_pain_mol()

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


def read_gene_list(infile, mode):
    gene_list = []
    try:
        with open(infile, mode) as FH:
            for eachGene in FH:
                gene_list.append(eachGene.strip("\n"))
        return gene_list

    except FileNotFoundError:
        sys.exit("Given file doesn't exist")

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
    # return set(filtered_markers)
    return filtered_dict


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
