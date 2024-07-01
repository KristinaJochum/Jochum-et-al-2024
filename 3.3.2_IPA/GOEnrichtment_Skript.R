##get genes of interest per substance in one file 
###genes of interest are: p < 0.05 (not adjusted) and Fold regulation < -1.5 or >1.5

#read excel files

library(readr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(plotly)
library(ggpubr)
library(cowplot)

Results_HepaRG <- as.data.frame(read_delim("HepaRG_IPA_Upload.csv", 
                                           delim = ";", escape_double = FALSE, 
                                           trim_ws = TRUE))

Results_RPTEC <- as.data.frame(read_delim("RPTEC_IPA_Upload.csv", 
                                          delim = ";", escape_double = FALSE, 
                                          trim_ws = TRUE, n_max = 84))

Results_HepaRG_RPTEC <- as.data.frame(read_delim("HepaRG_RPTEC_IPA_Upload.csv", 
                                                 delim = ";", escape_double = FALSE, 
                                                 trim_ws = TRUE, n_max = 454))

#genes of interest
# find interesting genes

interesting_genes_function <- function(df_results){
  for(i in seq(2, 13, 2)){
    interesting_genes <- mutate(df_results[, c(i, i+1)],
                                interesting_candidate_fc=ifelse(df_results[i] > -1.5 & df_results[i] < 1.5,"NO","YES"),
                                interesting_candidate_pval=ifelse(df_results[i+1] >= 0.05,"NO","YES"),
                                interesting_candidate=ifelse(interesting_candidate_fc=="YES"&interesting_candidate_pval=="YES","YES","NO"),
                                .keep = "none")
    names(interesting_genes) <- gsub("interesting_candidate", paste0("GO_genes_", gsub("_Fold_Regulation", "", names(df_results[i]))), names(interesting_genes))
    interesting_genes <- apply(interesting_genes, 2, unlist)
    df_results <- data.frame(df_results, as.data.frame(interesting_genes))
  }
  return(df_results)
}

interesting_genes_HepaRG <- interesting_genes_function(Results_HepaRG)
interesting_genes_RPTEC <- interesting_genes_function(Results_RPTEC)
interesting_genes_HepaRG_RPTEC <- interesting_genes_function(Results_HepaRG_RPTEC)


#r create GORilla lists
##input for function: df from above, CellLine in ""

GORilla_Lists <- function(interesting_genes_df, CellLine){
  for(c in seq(16, 31, 3)){
    if(sum(interesting_genes_df[,c]=="YES")>0){
      write_delim(interesting_genes_df %>% filter(interesting_genes_df[,c] == "YES") %>% select(Refseq),
                  file=paste0("GORilla_input_", CellLine,
                              gsub("GO_genes_", "_", names(interesting_genes_df[c])),
                              ".tsv"),
                  col_names = FALSE)
      write_delim(interesting_genes_df %>% filter(interesting_genes_df[,c] == "NO") %>% select(Refseq),
                  file = paste0("GOrilla_input_Symbol_Background_", CellLine,
                                gsub("GO_genes_", "_", names(interesting_genes_df[c])),
                                ".tsv"),
                  col_names = F)
      
      print(paste(gsub("GO_genes_", "", names(interesting_genes_df[c])),"list of characteristically changed genes written"))
      
    }else{
      print(paste(gsub("GO_genes_", "", names(interesting_genes_df[c])),"doesn't have characteristically changed genes"))
    }
  }
}

GORilla_Lists(interesting_genes_HepaRG, "HepaRG")
GORilla_Lists(interesting_genes_HepaRG_RPTEC, "HepaRG_RPTEC")
GORilla_Lists(interesting_genes_RPTEC, "RPTEC")

