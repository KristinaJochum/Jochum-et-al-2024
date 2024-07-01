#Array - data from GeneGlobe
#HepaRG

library(readr)

HepaRG_IPA_Upload <- as.data.frame(read_delim("HepaRG_IPA_Upload.csv", 
                                    delim = ";", escape_double = FALSE, 
                                    trim_ws = TRUE))

HepaRG_Array_Genes <- read_delim("HepaRG_Genes_Molecular_Toxicology_Pathway_Finder.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

row.names(HepaRG_IPA_Upload) <- make.names(HepaRG_Array_Genes$Symbol[1:370], unique = T)

library(reshape)

HepaRG_Array_genes_per_effect <- read_delim("HepaRG_genes_per_effect_MolToxPathFind.csv", 
                                            delim = ";", 
                                            escape_double = FALSE, 
                                            trim_ws = TRUE)

HepaRG_IPA_Upload$Gene <- gsub("HLA.DRB1", "HLA-DRB1", row.names(HepaRG_IPA_Upload))

HepaRG_IPA_Upload_effect <- merge(
  x = HepaRG_Array_genes_per_effect, 
  y = HepaRG_IPA_Upload, 
  by.x = "Gene", 
  by.y = "Gene", 
  all = T)

HepaRG_IPA_Upload_effect <- HepaRG_IPA_Upload_effect[order(
  HepaRG_IPA_Upload_effect$effect, 
  HepaRG_IPA_Upload_effect$Gene, 
  decreasing = c(F,F), 
  method = "radix"),]

row.names(HepaRG_IPA_Upload_effect) <- paste(
  HepaRG_IPA_Upload_effect$effect, 
  HepaRG_IPA_Upload_effect$Gene, 
  sep = "_")


##complexHeatmap

library(ComplexHeatmap)

###inlcude levels in dataframe to be used for complexHeatmap

for(i in grep("_Fold_Regulation", names(HepaRG_IPA_Upload_effect), ignore.case = T)){
  levels <- cut(HepaRG_IPA_Upload_effect[,i], 
                breaks = c(-Inf, -2, -1.5, 1.50, 2.00, Inf), 
                labels = c("< 0.50", "0.50-0.75", "0.75-1.50", "1.50-2.00", "> 2.00"),
                right = FALSE)
  HepaRG_IPA_Upload_effect <-  cbind(HepaRG_IPA_Upload_effect, levels)
  names(HepaRG_IPA_Upload_effect)[ncol(HepaRG_IPA_Upload_effect)] <- paste0(names(HepaRG_IPA_Upload_effect)[i], "_levels")
}

for(i in grep("_pvalue", names(HepaRG_IPA_Upload_effect), ignore.case = T)){
  asterisks <- cut(HepaRG_IPA_Upload_effect[,i], 
                breaks = c(0, 0.05, Inf), 
                labels = c("*", ""),
                right = FALSE)
  HepaRG_IPA_Upload_effect <-  cbind(HepaRG_IPA_Upload_effect, asterisks)
  names(HepaRG_IPA_Upload_effect)[ncol(HepaRG_IPA_Upload_effect)] <- paste0(names(HepaRG_IPA_Upload_effect)[i], "_asterisks")
}


### set colors for heatmap

colors <- structure(
  c("dodgerblue2", "lightskyblue1", "white", "lightpink", "red2"), 
  names = c("< 0.50", "0.50-0.75", "0.75-1.50", "1.50-2.00", "> 2.00"))

### adapt row and col titles 

HepaRG_IPA_Upload_effect$effect <-  gsub("%", "\n", HepaRG_IPA_Upload_effect$effect)
col_title <- c("Cyproconazole 40 µM","Fluxapyroxad 20 µM","Azoxystrobin 15 µM","Chlorotoluron 250 µM","Thiabendazole 100 µM","2-Phenylphenol 80 µM")

###Heatmap

HT_1 <- Heatmap(HepaRG_IPA_Upload_effect[c(1:137), c(grep("_levels", names(HepaRG_IPA_Upload_effect), ignore.case = T))],
                width = unit(20, "mm"),
                show_heatmap_legend = F,
                name = "Fold Change",
                col = colors,
                column_names_rot = 90,
                row_labels = HepaRG_IPA_Upload_effect$Gene[c(1:137)],
                row_names_gp = gpar(fontsize = 6),
                split = HepaRG_IPA_Upload_effect$effect[c(1:137)],
                row_title_rot = 360,
                row_title_gp = gpar(fontsize = 7, just = "center"),
                column_labels = col_title,
                column_names_centered = F,
                column_names_gp = gpar(fontsize = 7),
                cluster_columns = F, 
                cluster_rows = F, 
                row_names_side = "right",
                row_title_side = "right",
                border = T,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(HepaRG_IPA_Upload_effect[i, j+21], x, y, gp= gpar(fontsize = 6))
                })

HT_1

HT_2 <- Heatmap(HepaRG_IPA_Upload_effect[c(138:267), c(grep("_levels", names(HepaRG_IPA_Upload_effect), ignore.case = T))],
                width = unit(20, "mm"),
                name = "Fold Change",
                col = colors,
                column_names_rot = 90,
                row_labels = HepaRG_IPA_Upload_effect$Gene[c(138:267)],
                row_names_gp = gpar(fontsize = 6),
                row_title_rot = 360,
                row_title_gp = gpar(fontsize = 7, just = "center"),
                column_labels = col_title,
                column_names_centered = F,
                column_names_gp = gpar(fontsize = 7),
                cluster_columns = F, 
                cluster_rows = F, 
                row_names_side = "right", 
                row_title_side = "right",
                split = HepaRG_IPA_Upload_effect$effect[c(138:267)],
                border = T,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(HepaRG_IPA_Upload_effect[i+137,j+21], x, y, gp= gpar(fontsize = 6))
                },
                show_heatmap_legend = F)

HT_2

HT_3 <- Heatmap(HepaRG_IPA_Upload_effect[c(268:398), c(grep("_levels", names(HepaRG_IPA_Upload_effect), ignore.case = T))],
                width = unit(20, "mm"),
                name = "Fold Change",
                col = colors,
                column_names_rot = 90,
                row_labels = HepaRG_IPA_Upload_effect$Gene[c(268:398)],
                row_names_gp = gpar(fontsize = 6),
                row_title_rot = 360,
                row_title_gp = gpar(fontsize = 7, just = "center"),
                column_labels = col_title,
                column_names_centered = F,
                column_names_gp = gpar(fontsize = 7),
                cluster_columns = F, 
                cluster_rows = F, 
                row_names_side = "right", 
                row_title_side = "right",
                split = HepaRG_IPA_Upload_effect$effect[c(268:398)],
                border = T,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(HepaRG_IPA_Upload_effect[i+267,j+21], x, y, gp= gpar(fontsize = 6))
                })

HT_3

###grid of heatmaps

png(filename = "HepaRG_Array.png", 
    width = 220,
    height = 280,
    units = "mm",
    res = 900)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 3,
                                           widths = unit(c(0.9, 0.9, 1.2), "null"),
                                           heights = unit(1, "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(HT_1, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(HT_2, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(HT_3, newpage = FALSE)
upViewport()

dev.off()

##count times per effect when FR>1.5 or FC<-1.5 + significant

df <- HepaRG_IPA_Upload_effect[, order(names(HepaRG_IPA_Upload_effect))]
effects <- unique(HepaRG_IPA_Upload_effect$effect)
new_df <- data.frame(0)

for(k in grep("_Fold_regulation$", names(df), ignore.case = T)){
  affected <- numeric(0)
  for(m in 1:nrow(df)){
    if(df[m,k] > 1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
    if(df[m,k] < -1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
  }
  percent_vector <- numeric(0)
  for (i in effects) {
    Summe <- numeric(0)
    Summe <- length(grep(i, affected))
    total <- nrow(df[c(grep(i, df$effect)),])
    Percent <- round(Summe*100/total)
    percent_vector <- c(percent_vector, Percent)
  }
  new_df <- cbind(new_df, percent_vector)
}

names(new_df)[2:7] <- names(df[grep("_Fold_regulation$", names(df), ignore.case = T)])
row.names(new_df) <- effects
new_df$effects <- effects

##max value
HepaRG_IPA_Upload["CYP2B6", c(grep("Fold_Regulation", names(HepaRG_IPA_Upload), ignore.case = T))]
max_value <- max(HepaRG_IPA_Upload[c(grep("Fold_Regulation", names(HepaRG_IPA_Upload), ignore.case = T))])
max_position <- which(df == max_value, arr.ind = TRUE)
row_name <- rownames(df)[max_position[1]]
col_name <- colnames(df)[max_position[2]]
cat("The highest value is", max_value, "located at row", row_name, "and column", col_name, "\n")


##Supplement

HepaRG_Array_rep1 <- as.data.frame(read_delim("HepaRG_Array_rep1_data.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE))

HepaRG_Array_rep2 <- as.data.frame(read_delim("HepaRG_Array_rep2_data.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE))

HepaRG_Array_rep3 <- as.data.frame(read_delim("HepaRG_Array_rep3_data.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE))

cutoff <- function(dataframe_column){
  for(i in 1:length(dataframe_column)){
    if(is.na(dataframe_column[i]) == TRUE){
      dataframe_column[i] <- 35
    }
    if(dataframe_column[i] > 35){
      dataframe_column[i] <- 35
    }
  }
  return(dataframe_column[c(1:length(dataframe_column))])
}

HepaRG_Array_rep1 <- as.data.frame(apply(HepaRG_Array_rep1, 2, cutoff))
HepaRG_Array_rep2 <- as.data.frame(apply(HepaRG_Array_rep2, 2, cutoff))
HepaRG_Array_rep3 <- as.data.frame(apply(HepaRG_Array_rep3, 2, cutoff))

###sheet including gene names in same order as on array

HepaRG_Array_Genes <- read_delim("HepaRG_Genes_Molecular_Toxicology_Pathway_Finder.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

###assign gene names as rownames

row.names(HepaRG_Array_rep1) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)
row.names(HepaRG_Array_rep2) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)
row.names(HepaRG_Array_rep3) <- make.names(HepaRG_Array_Genes$Symbol, unique = T)


###function (applicable to all reps) for -ddct
###HKG mean NC
###dCt NC
###HKG mean Substances

calculate_ddCt <- function(df_Ct){
  HKG_mean_NC <- mean(df_Ct[c(371:375),1])
  dCt_NC <- numeric(0)
  for(k in 1:nrow(df_Ct)){
    dCt_NC_gene <- df_Ct[k,1] - HKG_mean_NC
    dCt_NC <- c(dCt_NC, dCt_NC_gene)
  }
  df_ddCt <- numeric(0)
  for(i in 2:ncol(df_Ct)){
    HKG_mean_substance <- mean(df_Ct[c(371:375),i])
    ddCt_substance <- numeric(0)
    for(m in 1:370){
      dCt_substance_gene <- df_Ct[m,i] - HKG_mean_substance
      ddCt_substance_gene <- dCt_substance_gene - dCt_NC[m]
      ddCt_substance <- c(ddCt_substance, ddCt_substance_gene)
    }
    df_ddCt <- as.data.frame(cbind(df_ddCt, ddCt_substance))
    colnames(df_ddCt)[i-1] <- colnames(df_Ct)[i]
  }
  row.names(df_ddCt) <- row.names(df_Ct)[1:370]
  return(df_ddCt)
}

HepaRG_Array_ddCt_rep1 <- calculate_ddCt(HepaRG_Array_rep1)
HepaRG_Array_ddCt_rep2 <- calculate_ddCt(HepaRG_Array_rep2)
HepaRG_Array_ddCt_rep3 <- calculate_ddCt(HepaRG_Array_rep3)

colnames(HepaRG_Array_ddCt_rep1) <- paste(colnames(HepaRG_Array_ddCt_rep1), "rep1", sep = "_")
colnames(HepaRG_Array_ddCt_rep2) <- paste(colnames(HepaRG_Array_ddCt_rep2), "rep2", sep = "_")
colnames(HepaRG_Array_ddCt_rep3) <- paste(colnames(HepaRG_Array_ddCt_rep3), "rep3", sep = "_")

HepaRG_Array_ddCt_allrep <- cbind(HepaRG_Array_ddCt_rep1, HepaRG_Array_ddCt_rep2, HepaRG_Array_ddCt_rep3)

library(ggplot2)

boxplot_Array_HepaRG <- function(Substance){
  names(HepaRG_Array_ddCt_allrep) <- gsub("X2.", "2-", names(HepaRG_Array_ddCt_allrep))
  Sub <- grep(Substance, names(HepaRG_Array_ddCt_allrep))
  datamelt_boxplot <- melt(cbind(Genes = row.names(HepaRG_Array_ddCt_allrep), HepaRG_Array_ddCt_allrep[, Sub]))
  datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", datamelt_boxplot$variable)
  print(unique(datamelt_boxplot$variable))
  datamelt_boxplot$value <- (-1)*datamelt_boxplot$value
  col_title <- c("Cyproconazole 40 µM","Fluxapyroxad 20 µM","Azoxystrobin 15 µM","Chlorotoluron 250 µM","Thiabendazole 100 µM","2-Phenylphenol 80 µM")
  datamelt_boxplot$variable <- gsub(Substance, col_title[grep(Substance, col_title)],datamelt_boxplot$variable)
  print(unique(datamelt_boxplot$variable))
  Plot <- ggplot(datamelt_boxplot, aes(Genes, value))+
    geom_point(size = 1)+
    geom_boxplot(width = 0.5, linewidth = 0.3)+
    stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.2,
      linewidth = 0.3
    )+
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 1
    )+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5)+
    geom_hline(yintercept = -1, linetype = "dashed", color = "red", linewidth = 0.5)+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 6),
      strip.background = element_blank(), 
      strip.text.x = element_blank())+
    ylab(expression("relative gene transcription "*log[2]*"FC"))+
    facet_wrap(~ ifelse(as.numeric(as.factor(Genes)) <= 185, "Oben", "Unten"), scales = "free", ncol = 1) +
    ggtitle(paste("HepaRG", "\u2012", unique(datamelt_boxplot$variable), " "))
  return(Plot)
}

Pesticides <- c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "2-Phenylphenol")

for(i in Pesticides){
  p <- boxplot_Array_HepaRG(i)
  ggsave(
    filename = paste0(i, ".png"),
    plot = p,
    device = "png",
    width = 297,
    height = 210,
    units = "mm",
    dpi = 900)
}



##RPTEC

library(readr)

RPTEC_IPA_Upload <- as.data.frame(read_delim("RPTEC_IPA_Upload.csv", 
                                              delim = ";", escape_double = FALSE, 
                                              trim_ws = TRUE, n_max = 84))

RPTEC_Array_Genes <- read_delim("RPTEC_Array_GeneList.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

row.names(RPTEC_IPA_Upload) <- make.names(RPTEC_Array_Genes$`Gene Symbol`[1:84], unique = T)


library(reshape)

RPTEC_Array_genes_per_effect <- read_delim("RPTEC_genes_per_effect_nephotoxicity.csv", 
                                            delim = ";", 
                                            escape_double = FALSE, 
                                            trim_ws = TRUE)

RPTEC_IPA_Upload$Gene <- row.names(RPTEC_IPA_Upload)

RPTEC_IPA_Upload_effect <- merge(
  x = RPTEC_Array_genes_per_effect, 
  y = RPTEC_IPA_Upload, 
  by.x = "Gene", 
  by.y = "Gene", 
  all = T)

RPTEC_IPA_Upload_effect <- RPTEC_IPA_Upload_effect[order(
  RPTEC_IPA_Upload_effect$effect, 
  RPTEC_IPA_Upload_effect$Gene, 
  decreasing = c(F,F), 
  method = "radix"),]

row.names(RPTEC_IPA_Upload_effect) <- paste(
  RPTEC_IPA_Upload_effect$effect, 
  RPTEC_IPA_Upload_effect$Gene, 
  sep = "_")


##complexHeatmap

library(ComplexHeatmap)

###inlcude levels in dataframe to be used for complexHeatmap

for(i in grep("_Fold_Regulation", names(RPTEC_IPA_Upload_effect), ignore.case = T)){
  levels <- cut(RPTEC_IPA_Upload_effect[,i], 
                breaks = c(-Inf, -2, -1.5, 1.50, 2.00, Inf), 
                labels = c("< 0.50", "0.50-0.75", "0.75-1.50", "1.50-2.00", "> 2.00"),
                right = FALSE)
  RPTEC_IPA_Upload_effect <-  cbind(RPTEC_IPA_Upload_effect, levels)
  names(RPTEC_IPA_Upload_effect)[ncol(RPTEC_IPA_Upload_effect)] <- paste0(names(RPTEC_IPA_Upload_effect)[i], "_levels")
}

for(i in grep("_pvalue", names(RPTEC_IPA_Upload_effect), ignore.case = T)){
  asterisks <- cut(RPTEC_IPA_Upload_effect[,i], 
                   breaks = c(0, 0.05, Inf), 
                   labels = c("*", ""),
                   right = FALSE)
  RPTEC_IPA_Upload_effect <-  cbind(RPTEC_IPA_Upload_effect, asterisks)
  names(RPTEC_IPA_Upload_effect)[ncol(RPTEC_IPA_Upload_effect)] <- paste0(names(RPTEC_IPA_Upload_effect)[i], "_asterisks")
}


### set colors for heatmap

colors <- structure(
  c("dodgerblue2", "lightskyblue1", "white", "lightpink", "red2"), 
  names = c("< 0.50", "0.50-0.75", "0.75-1.50", "1.50-2.00", "> 2.00"))

### adapt row and col titles 

RPTEC_IPA_Upload_effect$effect <-  gsub("%", "\n", RPTEC_IPA_Upload_effect$effect)
col_title <- c("Cyproconazole 300 µM","Fluxapyroxad 30 µM","Azoxystrobin 3 µM","Chlorotoluron 900 µM","Thiabendazole 300 µM","2-Phenylphenol 210 µM")

###Heatmap

png(filename = "RPTEC_Array.png", 
    width = 70,
    height = 280,
    units = "mm",
    res = 900)

Heatmap(RPTEC_IPA_Upload_effect[, c(grep("_levels", names(RPTEC_IPA_Upload_effect)))],
        name = "Fold Change",
        width = unit(20, "mm"),
        col = colors,
        column_names_rot = 90,
        row_labels = RPTEC_IPA_Upload_effect$Gene,
        row_names_gp = gpar(fontsize = 6),
        row_title = gsub(" ", "\n", unique(RPTEC_IPA_Upload_effect$effect)),
        row_title_rot = 360,
        row_title_gp = gpar(fontsize = 7, just = "center"),
        column_labels = col_title,
        column_names_centered = F,
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = F, 
        cluster_rows = F, 
        row_names_side = "right",
        row_title_side = "right",
        split = RPTEC_IPA_Upload_effect$effect,
        border = T,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize =7)),
        cell_fun = function(j, i, x, y, width, height, fill){
          grid.text(RPTEC_IPA_Upload_effect[i,j+21], x, y, gp= gpar(fontsize = 6))
        })

dev.off()

##count times per effect when FE>1.5 or FC<-1.5 + significant

df <- RPTEC_IPA_Upload_effect[, order(names(RPTEC_IPA_Upload_effect))]
effects <- unique(RPTEC_IPA_Upload_effect$effect)
new_df <- data.frame(0)

for(k in grep("_Fold_regulation$", names(df), ignore.case = T)){
  affected <- numeric(0)
  for(m in 1:nrow(df)){
    if(df[m,k] > 1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
    if(df[m,k] < -1.5 & df[m, k+3] == "*"){
      affected <- c(affected, df$effect[m])
    }
  }
  percent_vector <- numeric(0)
  for (i in effects) {
    Summe <- numeric(0)
    Summe <- length(grep(i, affected))
    total <- nrow(df[c(grep(i, df$effect)),])
    Percent <- round(Summe*100/total)
    percent_vector <- c(percent_vector, Percent)
  }
  new_df <- cbind(new_df, percent_vector)
}

names(new_df)[2:7] <- names(df[grep("_Fold_regulation$", names(df), ignore.case = T)])
row.names(new_df) <- effects
new_df$effects <- effects

##max value
RPTEC_IPA_Upload["IGFBP1", c(grep("Fold_Regulation", names(RPTEC_IPA_Upload), ignore.case = T))]
max_value <- max(RPTEC_IPA_Upload[c(grep("Fold_Regulation", names(RPTEC_IPA_Upload), ignore.case = T))])
max_position <- which(df == max_value, arr.ind = TRUE)
row_name <- rownames(df)[max_position[1]]
col_name <- colnames(df)[max_position[2]]
cat("The highest value is", max_value, "located at row", row_name, "and column", col_name, "\n")

##Supplementary

RPTEC_Array_rep1 <- as.data.frame(read_delim("RPTEC_Array_rep1_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

RPTEC_Array_rep2 <- as.data.frame(read_delim("RPTEC_Array_rep2_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

RPTEC_Array_rep3 <- as.data.frame(read_delim("RPTEC_Array_rep3_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))

###cutoff every Ct above 35

cutoff <- function(dataframe_column){
  for(i in 1:length(dataframe_column)){
    if(is.na(dataframe_column[i]) == TRUE){
      dataframe_column[i] <- 35
    }
    if(dataframe_column[i] > 35){
      dataframe_column[i] <- 35
    }
  }
  return(dataframe_column[c(1:length(dataframe_column))])
}

RPTEC_Array_rep1 <- as.data.frame(apply(RPTEC_Array_rep1, 2, cutoff))
RPTEC_Array_rep2 <- as.data.frame(apply(RPTEC_Array_rep2, 2, cutoff))
RPTEC_Array_rep3 <- as.data.frame(apply(RPTEC_Array_rep3, 2, cutoff))

###sheet including gene names in same order as on array

RPTEC_Array_GeneList <- read_csv("D:/BfR/Labor/Results/RPTEC/Array/.txt_data/Pesticides/RPTEC_Array_GeneList.csv")

###assign gene names as rownames

row.names(RPTEC_Array_rep1) <- make.names(RPTEC_Array_GeneList$`Gene Symbol`, unique = T)
row.names(RPTEC_Array_rep2) <- make.names(RPTEC_Array_GeneList$`Gene Symbol`, unique = T)
row.names(RPTEC_Array_rep3) <- make.names(RPTEC_Array_GeneList$`Gene Symbol`, unique = T)

calculate_ddCt <- function(df_Ct){
  HKG_mean_NC <- mean(df_Ct[c((nrow(df_Ct) - 11):(nrow(df_Ct) - 7)),1])
  dCt_NC <- numeric(0)
  for(k in 1:nrow(df_Ct)){
    dCt_NC_gene <- df_Ct[k,1] - HKG_mean_NC
    dCt_NC <- c(dCt_NC, dCt_NC_gene)
  }
  df_ddCt <- numeric(0)
  for(i in 2:ncol(df_Ct)){
    HKG_mean_substance <- mean(df_Ct[c((nrow(df_Ct) - 11):(nrow(df_Ct) - 7)),i])
    ddCt_substance <- numeric(0)
    for(m in 1:(nrow(df_Ct) - 12)){
      dCt_substance_gene <- df_Ct[m,i] - HKG_mean_substance
      ddCt_substance_gene <- dCt_substance_gene - dCt_NC[m]
      ddCt_substance <- c(ddCt_substance, ddCt_substance_gene)
    }
    df_ddCt <- as.data.frame(cbind(df_ddCt, ddCt_substance))
    colnames(df_ddCt)[i-1] <- colnames(df_Ct)[i]
  }
  row.names(df_ddCt) <- row.names(df_Ct)[1:(nrow(df_Ct) - 12)]
  return(df_ddCt)
}

RPTEC_Array_ddCt_rep1 <- calculate_ddCt(RPTEC_Array_rep1)
RPTEC_Array_ddCt_rep2 <- calculate_ddCt(RPTEC_Array_rep2)
RPTEC_Array_ddCt_rep3 <- calculate_ddCt(RPTEC_Array_rep3)

colnames(RPTEC_Array_ddCt_rep1) <- paste(colnames(RPTEC_Array_ddCt_rep1), "rep1", sep = "_")
colnames(RPTEC_Array_ddCt_rep2) <- paste(colnames(RPTEC_Array_ddCt_rep2), "rep2", sep = "_")
colnames(RPTEC_Array_ddCt_rep3) <- paste(colnames(RPTEC_Array_ddCt_rep3), "rep3", sep = "_")

RPTEC_Array_ddCt_allrep <- cbind(RPTEC_Array_ddCt_rep1, RPTEC_Array_ddCt_rep2, RPTEC_Array_ddCt_rep3)

library(ggplot2)

boxplot_Array_RPTEC <- function(pesticide1, pesticide2, save){
  names(RPTEC_Array_ddCt_allrep) <- gsub("X2.", "2-", names(RPTEC_Array_ddCt_allrep))
  Pesticides1_2 <- c(grep(pesticide1, names(RPTEC_Array_ddCt_allrep)), grep(pesticide2, names(RPTEC_Array_ddCt_allrep)))
  datamelt_boxplot <- melt(cbind(Genes = row.names(RPTEC_Array_ddCt_allrep), RPTEC_Array_ddCt_allrep[, Pesticides1_2]))
  datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", datamelt_boxplot$variable)
  print(unique(datamelt_boxplot$variable))
  datamelt_boxplot$value <- (-1)*datamelt_boxplot$value
  col_title <- c("Cyproconazole 300 µM","Fluxapyroxad 30 µM","Azoxystrobin 3 µM","Chlorotoluron 900 µM","Thiabendazole 300 µM","2-Phenylphenol 210 µM")
  datamelt_boxplot$variable <- ifelse(datamelt_boxplot$variable == pesticide1, col_title[grep(pesticide1, col_title)], col_title[grep(pesticide2, col_title)])
  print(unique(datamelt_boxplot$variable))
  Plot <- ggplot(datamelt_boxplot, aes(Genes, value))+
    geom_point(size = 1)+
    geom_boxplot(width = 0.5, linewidth = 0.3)+
    stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.2,
      linewidth = 0.3
    )+
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 1
    )+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5)+
    geom_hline(yintercept = -1, linetype = "dashed", color = "red", linewidth = 0.5)+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 8))+
    ylab(expression("relative gene transcription "*log[2]*"FC"))+
    facet_wrap(~variable, 
               ncol = 1)+
    ggtitle("RPTEC")
  if(save == "YES"){
    print("save")
    print(getwd())
    ggsave(
      filename = paste0(pesticide1, "_", pesticide2, ".png"),
      plot = Plot,
      device = "png",
      width = 297,
      height = 210,
      units = "mm",
      dpi = 900)
  }
  return(Plot)
}

boxplot_Array_RPTEC("Cyproconazole", "Fluxapyroxad", "YES")

boxplot_Array_RPTEC("Azoxystrobin", "2-Phenylphenol", "YES")

boxplot_Array_RPTEC("Chlorotoluron", "Thiabendazole", "YES")
