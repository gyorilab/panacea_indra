# Libraries
library(tm)
library(proustr)
library(tidyverse)
library(hrbrthemes)
library(VennDiagram)
library(RColorBrewer)

# Functions
# Chart
draw_venn <- function(set_list, file_loc, title){
  
  venn.diagram(
    x = set_list,
    category.names = c("Human Primary" , "iPSC " , "Mouse Primary"),
    filename = paste0(file_loc),
    output=T,
    main = title,
    
    # Output features
    imagetype="tiff" ,
    height = 1000 , 
    width = 1000 , 
    resolution = 250,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.5,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(0, 1, 2),
    #cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    #rotation = 1
  )
  
  
}

# Read the datasets
setwd("~/PycharmProjects/INDRA/Panacea/")
ALL_HP <- readxl::read_xlsx("Data/ALL HP (human primary) selective.xlsx")
ALL_HP <- ALL_HP[!is.na(ALL_HP$`Type(s)`), ]
ALL_ID <- readxl::read_xlsx("Data/ALL ID (iPSC derived) selective.xlsx")
ALL_ID <- ALL_ID[!is.na(ALL_ID$`Type(s)`), ]
ALL_MP <- readxl::read_xlsx("Data/ALL MP (mouse primary) selective.xlsx")
ALL_MP <- ALL_MP[!is.na(ALL_MP$`Type(s)`), ]
# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")

# Creating generalized sets
set_ALL_HP <- na.omit(ALL_HP$Symbol)
set_ALL_ID <- na.omit(ALL_ID$Symbol)
set_ALL_MP <- na.omit(ALL_MP$Symbol)


all_datasets <- list(
  'Human Primary' = ALL_HP,
  'iPSC' = ALL_ID,
  'Mouse Primary' = ALL_MP
)

# Take a list of all the available cell types
all_types <- names(table(ALL_HP$`Type(s)`))

for(type in all_types){
  count = 1
  all_markers <- list()
  for(dataset in all_datasets){
    markers = c()
    for(s in 1:nrow(dataset)){
      if(dataset[s,5]==type){
        markers <- c(markers, dataset$Symbol[s])
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
