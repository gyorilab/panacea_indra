library(tximeta)
library(stringr)
library(tximport)


wd <- './gitHub/panacea_indra/pain_model/data/Primary_mouse/salmon_GRCm38.p6/'
setwd(wd)

# Make txdb
txdb <- makeTxDbFromGFF(file = '../gencode.vM24.primary_assembly.annotation.gtf')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files <- file.path(list.files('./'), "quant.sf")
s <- "(\\w+_TRA\\w+\\d+_)([\\w-\\d\\_]+.S\\d+)"

coldata <- data.frame(files, names=str_match(files, s)[,3], 
                      stringsAsFactors=FALSE)

coldata <- data.frame('names'=coldata)
all(file.exists(coldata$names.files))

files <- coldata$names.files
names(files) <- coldata$names.names

txi <- tximport(files, type = "salmon", 
                tx2gene = tx2gene, 
                ignoreAfterBar = T)
txi.count <- txi$counts
dir.create('../output/', showWarnings = F)
write.csv(txi.count, '../output/primary_mouse_tpm.csv')
