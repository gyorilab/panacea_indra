# Libraries
library(tm)
library(proustr)
library(tidyverse)
library(hrbrthemes)
library(VennDiagram)
library(RColorBrewer)

# Functions
# Chart
draw_venn <- function(set_list, file_loc, title, cat_names, mycol, cat_dist){
  
  venn.diagram(
    x = set_list,
    category.names = cat_names,
    filename = paste0(file_loc),
    output=T,
    main = title,
    
    # Output features
    imagetype="tiff" ,
    height = 1200, 
    width = 1500 , 
    resolution = 250,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = mycol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.5,
    cat.fontface = "bold",
    #cat.default.pos = "inner",
    cat.pos = c(0, 1),
    cat.dist = cat_dist,
    cat.fontfamily = "sans",
    #rotation = 1
  )
  
  
}

# Read the datasets
setwd("~/gitHub/panacea_indra/other/draw_venn/")
list_1 <- readxl::read_xlsx("Data/list_060121/List1 CommongenesRayDRGvsRayBrainSCHeart.xlsx",
                            col_names = F)
list_1 <- list_1[!is.na(list_1$...8), ]

list_2 <- readxl::read_xls("Data/list_060121/List2 CommongenesOurhDRGvsRayBrainSCHeart.xls",
                            col_names = F)
list_2 <- list_2[!is.na(list_2$...4), ]

list_3 <- readxl::read_xlsx("Data/list_060121/List3 Commongenes(1stexpt)mDRGvs(2ndexpt)mDRG,BrainSCHeart.xlsx",
                            col_names = F)
list_3 <- list_3[!is.na(list_3$...9), ]

list_4 <- readxl::read_xlsx("Data/list_060121/List4 CommongenesNoci4wvsCorticalMotorCardio.xlsx",
                            col_names = F)
list_4 <- list_4[!is.na(list_4$...7), ]

list_5 <- readxl::read_xlsx("Data/list_060121/List5 CommongenesNoci8wvsCorticalMotorCardio.xlsx",
                            col_names = F)
list_5 <- list_5[!is.na(list_5$...7), ]


# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(5, "Pastel2")

# Creating generalized sets
set_list_1 <- list_1[,c(1,8)]
set_list_2 <- list_2[,c(1,4)]
set_list_3 <- list_3[,c(1,9)]
set_list_4 <- list_4[,c(1,7)]
set_list_5 <- list_5[,c(1,7)]



# Comparison 1 (1 vs 3 vs 4)
type <- 'all_genes'
cat_dist <- rep(0.052,5)
dir.create(paste0('./venn_diagrams/', type),
           showWarnings = F, 
           recursive = T)

loc = paste0('./venn_diagrams/', type, '/1vs3vs4.tiff')
all_genes <- list(set_list_1$...1, set_list_3$...1, set_list_4$...1)
cat_names <- c('CommongenesRayDRGvsRayBrainSCHeart',
               'Commongenes(1stexpt)mDRGvs(2ndexpt)mDRG,BrainSCHeart',
               'CommongenesNoci4wvsCorticalMotorCardio'
               )
draw_venn(all_genes, loc, type, cat_names, myCol[1:3], cat_dist[1:3])

# Comparison 2 (1 vs 3 vs 5)
loc = paste0('./venn_diagrams/', type, '/1vs3vs5.tiff')
all_genes <- list(set_list_1$...1, set_list_3$...1, set_list_5$...1)
cat_names <- c('CommongenesRayDRGvsRayBrainSCHeart',
               'Commongenes(1stexpt)mDRGvs(2ndexpt)mDRG,BrainSCHeart',
               'CommongenesNoci8wvsCorticalMotorCardio'
)
draw_venn(all_genes, loc, type, cat_names, myCol[1:3], cat_dist[1:3])



# Comparison 2 (2 vs 3 vs 4)
loc = paste0('./venn_diagrams/', type, '/2vs3vs4.tiff')
all_genes <- list(set_list_2$...1, set_list_3$...1, set_list_4$...1)
cat_names <- c('CommongenesOurhDRGvsRayBrainSCHeart',
               'Commongenes(1stexpt)mDRGvs(2ndexpt)mDRG,BrainSCHeart',
               'CommongenesNoci4wvsCorticalMotorCardio')
draw_venn(all_genes, loc, type, cat_names, myCol[1:3], cat_dist[1:3])


# Comparison 2 (2 vs 3 vs 5)
loc = paste0('./venn_diagrams/', type, '/2vs3vs5.tiff')
all_genes <- list(set_list_2$...1, set_list_3$...1, set_list_5$...1)
cat_names <- c('CommongenesOurhDRGvsRayBrainSCHeart',
               'Commongenes(1stexpt)mDRGvs(2ndexpt)mDRG,BrainSCHeart',
               'CommongenesNoci8wvsCorticalMotorCardio')
draw_venn(all_genes, loc, type, cat_names, myCol[1:3], cat_dist[1:3])


# Comparison 2 (1 vs 2)
loc = paste0('./venn_diagrams/', type, '/1vs2.tiff')
all_genes <- list(set_list_1$...1, set_list_2$...1)
cat_names <- c('CommongenesRayDRGvsRayBrainSCHeart',
               'CommongenesOurhDRGvsRayBrainSCHeart')

draw_venn(all_genes, loc, type, cat_names, myCol[1:2], cat_dist[1:2])


# Comparison 2 (4 vs 5)
loc = paste0('./venn_diagrams/', type, '/4vs5.tiff')
all_genes <- list(set_list_4$...1, set_list_5$...1)
cat_names <- c('CommongenesNoci4wvsCorticalMotorCardio',
               'CommongenesNoci8wvsCorticalMotorCardio')

draw_venn(all_genes, loc, type, cat_names, myCol[1:2], cat_dist[1:2])

all_datasets <- list(
  'CommongenesRayDRGvsRayBrainSCHeart' = set_list_1,
  'CommongenesOurhDRGvsRayBrainSCHeart' = set_list_2,
  'Commongenes(1stexpt)mDRGvs(2ndexpt)mDRG,BrainSCHeart' = set_list_3,
  'CommongenesNoci4wvsCorticalMotorCardio' = set_list_4,
  'CommongenesNoci8wvsCorticalMotorCardio' = set_list_5
)

# Take a list of all the available cell types
#all_types <- Reduce(intersect, all_datasets)
al_types <- c("other", "transcription regulator", "transporter", "enzyme",
              "transmembrane receptor", "ion channel", 
              "G-protein coupled receptor", "kinase", "phosphatase" )

for(type in all_types){
  count = 1
  all_markers <- list()
  for(dataset in all_datasets){
    markers = c()
    for(s in 1:nrow(dataset)){
      if(dataset[s, 2]==type){
        markers <- c(markers, dataset[s,1])
      }
    }

    if(length(markers) > 0){
      all_markers[names(all_datasets)[count]][[1]] <- markers
    }else{
      all_markers[names(all_datasets)[count]][[1]] <- ""
    }
    count = count + 1
  }
  dir.create(paste0('./venn_diagrams/', type),
             showWarnings = F, 
             recursive = T)
  loc = paste0('./venn_diagrams/', type, '/venn_diagram.tiff')
  draw_venn(all_markers, loc, type)
}

#draw_venn(list(set_ALL_HP, set_ALL_ID, set_ALL_MP), file_loc = './venn_diagram.tiff', 'All markers')
