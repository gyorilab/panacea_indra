#!usr/bin/python
import csv
import sys
import os



#### variable ####
code_dir_name = "/Users/sbunga/gitHub/singlecell_model_analysis/"
directory = "/Users/sbunga/gitHub/singlecell_model_analysis/output/pseudo_bulk_analysis/diff_files"
output_directory = code_dir_name+"/output/pseudo_bulk_analysis/divenn_files/"

# where is this file???
gene2infoFile = code_dir_name+"/input_files/biomart_geneInfo_mus.txt"
gene2info = {}

############## read gene assocaited information from biomart table 
with open(gene2infoFile) as infile:
    next(infile)
    for line in infile:
        data = line.strip().split("\t")
        #print(len(data))
        #print line.strip()
        if not data[0] in gene2info:
            #print(data[0])
            if len(data)>2:
                #print(data[0])
                gene2info[data[0]] = data[1]+"\t"+data[2]
            else:
                gene2info[data[0]] = data[1]+"\t-"


gene2FCs={}
gene2FC_noversions={}
DEs={}
deseq_files=os.listdir(directory)
for deseq_file in deseq_files:
    print(deseq_file)
    filename_base = deseq_file.split(".")[0]
    diff_filename = output_directory + "/diff_" + filename_base + ".txt"
    divenn_filename = output_directory + "/divenn_" + filename_base + ".txt"
    pathview_filename = output_directory + "/pathview_" + filename_base + ".txt"
    gene2FC = {}
    gene2FC_noversion = {}
    DE = []
    with open(directory + "/" + deseq_file) as csvfile, open(diff_filename,"w") as difffile, open(divenn_filename, "w") as outfile, open(pathview_filename,"w") as pathviewfile:
        infile = csv.reader(csvfile, delimiter=',')

        head = next(infile)
        print(head)
        difffile.write("\t".join(head)+"\tSymbol\tDescription\n")
        pathviewfile.write("GeneID\t"+filename_base+"\n")
        for line in infile:
            data = line
            
            pathviewfile.write(data[0].split(".")[0]+"\t"+data[2]+"\n")
            if data[5] != "NA":
                if not data[0].split(".")[0] in gene2FC_noversion:
                    gene2FC_noversion[data[0].split(".")[0]] = data[3]
                if float(data[5]) < float(0.05):
                    if not data[0] in DE:
                        DE.append(data[0])
                    if not data[0] in gene2FC:
                        gene2FC[data[0]] = data[2]
                    outfile.write(data[0]+"\t")
                    if float(data[2]) >float(1):
                        outfile.write("1\n")
                    else:
                        outfile.write("2\n")
                    geneID = data[0].split(".")[0]
                    if geneID in gene2info:
                        difffile.write("\t".join(line) + "\t" + gene2info[geneID] + "\n")
                    else:
                        difffile.write("\t".join(line) + "\t-\t-\n")
    gene2FCs[filename_base]=gene2FC
    gene2FC_noversions[filename_base]=gene2FC_noversion
    DEs[filename_base]=DE
        
############## identify the common/overlapping genes between experimental comparisons 

DE_list=list(DEs.values())

# get the intersection

commonGene = set(DE_list[0]).intersection(*DE_list)
commonGeneFile = output_directory+"/Overlapping_DE_gene.txt"

with open(commonGeneFile,"w") as outfile:
    header=""
    conditions=gene2FCs.keys()
    for condition in conditions:
        header+="log2FoldChange_"+condition+"\t"
    outfile.write("GeneID\tSymbol\tDescription\t"+header+"\n")
    for gene in commonGene:
        geneID = gene.split(".")[0]
        values=""
        for condition in conditions:
            values+=gene2FCs[condition][gene]+"\t"
        if geneID in gene2info:
            outfile.write(gene+"\t" + gene2info[geneID]+"\t" + values +"\n")
        else:
            outfile.write(gene+"\t-\t-t" + values +"\n")

