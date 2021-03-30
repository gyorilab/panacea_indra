#!/usr/bin/env nextflow

params.input = "/Users/sbunga/PycharmProjects/INDRA/ligandReceptorInteractome/input/"
params.output = "/Users/sbunga/PycharmProjects/nextflow/output/"
params.pkl = "/Users/sbunga/PycharmProjects/nextflow/output/pkl"


goa_human = Channel.fromPath( "$params.input/goa_human.gaf" )

neuro_immune_data = Channel.fromPath( "$params.input/Neuroimmune gene list .xlsx" )

drugbank_pkl = Channel.fromPath( "$params.input/drugbank_5.1.pkl" )

ion_channels = Channel.fromPath( "$params.input/ion_channels.txt" )

surface_proteins = Channel.fromPath( "$params.input/Surface Proteins.xlsx" )

Human_Pain_Genes_DB = Channel.fromPath( "$params.input/Human_Pain_Genes_DB.tsv" )



process main_inputs {
    
    cache 'lenient'

    input:
  
    file GO_ANNOTATIONS from goa_human
    file DATA_SPREADSHEET from neuro_immune_data
    file DRUG_BANK_PKL from drugbank_pkl
    file ION_CHANNELS from ion_channels
    file SURFACE_PROTEINS_WB from surface_proteins
    file HUMAN_PAIN_DB from Human_Pain_Genes_DB
    

    output:
    file "targets_by_drug.pkl" into targets_by_drug
    file "all_enzymes.pkl" into all_enzymes
    file "full_ligand_set.pkl" into full_ligand_set


    script:
    """
    python3 $workflow.projectDir/scripts/process_main_inputs.py $params.input $params.output\
    targets_by_drug.pkl receptors_in_data.pkl all_enzymes.pkl full_ligand_set.pkl \
    $GO_ANNOTATIONS $DATA_SPREADSHEET $DRUG_BANK_PKL $ION_CHANNELS $SURFACE_PROTEINS_WB \
    $HUMAN_PAIN_DB
    """
}




process cell_types{
    
    cache 'lenient'

    input:
    file ALL_ENZYMES from all_enzymes
    file FULL_LIGAND_SET from full_ligand_set
    file TARGETS_BY_DRUG from targets_by_drug

    output:
    file '*_enzyme_product_dict.pkl' into enzyme_product_dict
    file '*_de_enzyme_product_interaction.pkl' into de_enzyme_product_interaction
    file '*_pain_interactions.pkl' into pain_interactions
    file '*_hashes_by_gene_pair.pkl' into hashes_by_gene_pair
    file '*_ligands_in_data.pkl' into ligands_in_data

    script:
    """
    python3 $workflow.projectDir/scripts/process_cell_types.py $params.input $params.output $ALL_ENZYMES \
    $FULL_LIGAND_SET $TARGETS_BY_DRUG
    """
}



process get_cell_type_indra_statements{
    
    cache 'lenient'

    input:
    file 'HASHES_BY_GENE_PAIR' from hashes_by_gene_pair
    file 'LIGANDS_IN_DATA' from ligands_in_data

    output:
    stdout result

    script:
    """
    python3 $workflow.projectDir/scripts/get_indra_statements.py --input $params.input --output $params.output --hashes_by_gene_pair HASHES_BY_GENE_PAIR* --ligands_in_data LIGANDS_IN_DATA*
    """

}

result.view({it.trim()})

/*
process test{
    
    cache 'lenient'

    input:
    file 'z' from enzyme_product_dict

    output:
    stdout result

    script:
    """
    python $workflow.projectDir/scripts/test.py z*
    
    """
    //python $workflow.projectDir/scripts/test.py ${x}

}
result.view({it.trim()})
*/



/*
process part_2 {

    cache 'lenient'

    input:
    file ALL_ENZYMES from all_enzymes
    file FULL_LIGAND_SET from full_ligand_set
    file RECEPTORS_IN_DATA from receptors_by_drug
    file TARGETS_BY_DRUG from targets_by_drug

    
    output:
    file "possible_en_drug_targets.pkl" into possible_en_drug_targets
    file "ligands_FC.pkl" into ligands_FC
    file "enzymes_FC.pkl" into enzymes_FC
    file "stmts_db_by_cell_type.pkl" into stmts_db_by_cell_type
    file "stmts_by_cell_type.pkl" into stmts_by_cell_type
    file "enzyme_product_dict.pkl" into enzyme_product_dict
    file "de_enzyme_products.pkl" into de_enzyme_products
    file "all_ligand_receptor_statements.pkl" into all_ligand_receptor_statements
    file "ligands_by_receptor.pkl" into ligands_by_receptor
    

    script:
    """
    python3 $workflow.projectDir/scripts/p2.py $params.input $params.output $ALL_ENZYMES \
    $FULL_LIGAND_SET $RECEPTORS_IN_DATA $TARGETS_BY_DRUG
    """
}
*/