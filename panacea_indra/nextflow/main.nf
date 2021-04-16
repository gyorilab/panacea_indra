#!/usr/bin/env nextflow

params.input = "/Users/sbunga/PycharmProjects/INDRA/ligandReceptorInteractome/input/"
params.output = "/Users/sbunga/gitHub/panacea_indra/panacea_indra/nextflow/output/"
params.pkl = "/Users/sbunga/gitHub/panacea_indra/panacea_indra/nextflow/output/pkl"


goa_human = Channel.fromPath( "$params.input/goa_human.gaf" )

neuro_immune_data = Channel.fromPath( "$params.input/Neuroimmune gene list .xlsx" )

drugbank_pkl = Channel.fromPath( "$params.input/drugbank_5.1.pkl" )

ion_channels = Channel.fromPath( "$params.input/ion_channels.txt" )

surface_proteins = Channel.fromPath( "$params.input/Surface Proteins.xlsx" )

Human_Pain_Genes_DB = Channel.fromPath( "$params.input/Human_Pain_Genes_DB.tsv" )

ligand_receptor_spreadsheet = Channel.fromPath( "$params.input/ncomms8866-s3.xlsx")


process main_inputs {
    
    cache 'lenient'

    input:
  
    file GO_ANNOTATIONS from goa_human
    file DATA_SPREADSHEET from neuro_immune_data
    file DRUG_BANK_PKL from drugbank_pkl
    file ION_CHANNELS from ion_channels
    file SURFACE_PROTEINS_WB from surface_proteins
    file LIGAND_RECEPTOR_SPREADSHEET from ligand_receptor_spreadsheet
    

    output:
    file "targets_by_drug.pkl" into targets_by_drug
    file "all_enzymes.pkl" into all_enzymes
    file "full_ligand_set.pkl" into full_ligand_set


    script:
    """
    python3 $workflow.projectDir/scripts/process_main_inputs.py $params.input $params.output\
    targets_by_drug.pkl receptors_in_data.pkl all_enzymes.pkl full_ligand_set.pkl \
    $GO_ANNOTATIONS $DATA_SPREADSHEET $DRUG_BANK_PKL $ION_CHANNELS $SURFACE_PROTEINS_WB \
    $LIGAND_RECEPTOR_SPREADSHEET
    """
}


targets_by_drug.into {td1; td2; td3; td4}


process cell_types{
    
    cache 'lenient'

    input:
    file ALL_ENZYMES from all_enzymes
    file FULL_LIGAND_SET from full_ligand_set
    file TARGETS_BY_DRUG from td1

    output:
    file 'enzyme_product_dict.pkl' into enzyme_product_dict
    file '*_de_enzyme_product_interaction.pkl' into de_enzyme_product_interaction
    file '*_pain_interactions.pkl' into pain_interactions
    file '*_hashes_by_gene_pair.pkl' into hashes_by_gene_pair
    file '*_ligands_in_data.pkl' into ligands_in_data
    file 'de_enzyme_product_list.pkl' into de_enzyme_product_list
    file 'possible_en_drug_targets.pkl' into possible_en_drug_targets
    file 'enzyme_possible_drug_targets.pkl' into enzyme_possible_drug_targets
    file 'ligands_FC.pkl' into ligands_FC
    file 'enzymes_FC.pkl' into enzymes_FC

    script:
    """
    python3 $workflow.projectDir/scripts/process_cell_types.py $params.input $params.output $ALL_ENZYMES \
    $FULL_LIGAND_SET $TARGETS_BY_DRUG
    """
}



ligands_in_data.into { lg_1; lg_2; lg_3; lg_4 }
enzyme_product_dict.into {en_pd_1; en_pd_2; en_pd_3}

process get_cell_type_indra_statements {
    
    cache 'lenient'

    input:
    file 'HASHES_BY_GENE_PAIR' from hashes_by_gene_pair
    file 'LIGANDS_IN_DATA' from lg_1
    file 'ENZYME_POSSIBLE_DRUG_TARGETS' from enzyme_possible_drug_targets

    output:
    file 'possible_drug_targets.pkl' into possible_drug_targets
    file 'possible_db_drug_targets.pkl' into possible_db_drug_targets
    file 'stmts_db_by_cell_type.pkl' into stmts_db_by_cell_type
    file 'stmts_by_cell_type.pkl' into stmts_by_cell_type
    file 'all_ligand_receptor_statements.pkl' into all_ligand_receptor_statements
    file '*_ligands_by_receptor.pkl' into ligands_by_receptor
    file '*_ligands_by_receptor_db.pkl' into ligands_by_receptor_db

    script:
    """
    python3 $workflow.projectDir/scripts/get_indra_statements.py --input $params.input --output $params.output --hashes_by_gene_pair HASHES_BY_GENE_PAIR* --ligands_in_data LIGANDS_IN_DATA* --enzyme_possible_drug_targets ENZYME_POSSIBLE_DRUG_TARGETS
    """

}



stmts_db_by_cell_type.into {stmt_db_cell_type_1; stmt_db_cell_type_2; stmt_db_cell_type_3;}


process plot_interaction_potential{
    cache 'lenient'

    input:
    file 'LIGANDS_IN_DATA' from lg_2
    file 'STMTS_DB_BY_CELL_TYPE' from stmt_db_cell_type_1

    output:


    script:
    """
    python3 $workflow.projectDir/scripts/plot_interaction_potential.py --input $params.input --output $params.output \
    --stmts_db_by_cell_type STMTS_DB_BY_CELL_TYPE --ligands_in_data LIGANDS_IN_DATA*
    """
}


possible_drug_targets.into {pdt1; pdt2}

process get_small_mol_report{
    
    cache 'lenient'

    input:
    file TARGETS_BY_DRUG from td2
    file 'POSSIBLE_DRUG_TARGETS' from pdt1
    file 'POSSIBLE_DB_DRUG_TARGETS' from possible_db_drug_targets

    output:

    script:
    """
    python3 $workflow.projectDir/scripts/get_small_mol_report.py $params.input $params.output  $TARGETS_BY_DRUG \
    $POSSIBLE_DRUG_TARGETS $POSSIBLE_DB_DRUG_TARGETS 
    """
}



process get_enzyme_interactions{
    
    cache 'lenient'

    input:
    file 'ALL_LIGAND_RECEPTOR_STATEMENTS' from all_ligand_receptor_statements
    file 'DE_ENZYME_PRODUCT_LIST' from de_enzyme_product_list

    output:
    file 'products_receptors.pkl' into products_receptors

    script:
    """
    python3 $workflow.projectDir/scripts/get_enzyme_interactions.py $params.input $params.output $DE_ENZYME_PRODUCT_LIST \
    $ALL_LIGAND_RECEPTOR_STATEMENTS
    """
}



process get_enzyme_product_drug_interactions{
    
    cache 'lenient'

    input:
    file POSSIBLE_EN_DRUG_TARGETS from possible_en_drug_targets
    file TARGETS_BY_DRUG from td3
    file ENZYME_PRODUCT_DICT from en_pd_1
    file 'LIGANDS_BY_RECEPTOR' from ligands_by_receptor


    output:


    script:
    """
    python3 $workflow.projectDir/scripts/get_enzyme_product_drug_interactions.py --input $params.input \
    --output $params.output \
    --possible_en_drug_targets $POSSIBLE_EN_DRUG_TARGETS \
    --targets_by_drug $TARGETS_BY_DRUG \
    --enzyme_product_dict $ENZYME_PRODUCT_DICT \
    --ligands_by_receptor LIGANDS_BY_RECEPTOR*
    """
}


process create_digraph{
    
    cache 'lenient'

    input:
    file ENZYME_PRODUCT_DICT from en_pd_2
    file PRODUCTS_RECEPTORS from products_receptors
    file STMT_DB_CELL_TYPE from stmt_db_cell_type_2
    file LIGANDS_FC from ligands_FC
    file ENZYMES_FC from enzymes_FC


    output:
    stdout result


    script:
    """
    python3 $workflow.projectDir/scripts/create_digraph.py --input $params.input \
    --output $params.output \
    --enzyme_product_dict $ENZYME_PRODUCT_DICT \
    --products_receptors $PRODUCTS_RECEPTORS \
    --stmts_db_by_cell_type $STMT_DB_CELL_TYPE \
    --enzymes_FC $ENZYMES_FC \
    --ligands_FC $LIGANDS_FC
    """
}

/*
process downstream_analysis{
    
    cache 'lenient'

    input:
    file HUMAN_PAIN_GENES_DB from Human_Pain_Genes_DB



    output:


    script:
    """
    python3 $workflow.projectDir/scripts/downstream_analysis.py --input $params.input \
    --output $params.output \
    --human_pain_db HUMAN_PAIN_GENES_DB
    """


}
*/

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



