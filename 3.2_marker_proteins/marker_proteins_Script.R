# Skript to evaluate Luminex data
##Import data
###import data

library(readr)
HepaRG_Protein_Result <- read_delim("HepaRG_Protein_Results.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE, col_select = c(1:15))

RPTEC_Protein_Result <- read_delim("RPTEC_Protein_Results.csv", 
                                                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

# create dataframe with rep1-3 side by side
##concentration as character for row.name
HepaRG_Protein_Result$concentration_uM <- as.character(HepaRG_Protein_Result$concentration_uM)

RPTEC_Protein_Result$concentration_uM <- as.character(RPTEC_Protein_Result$concentration_uM)

## one df per rep

HepaRG_Protein_Result_rep1 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 1),]
HepaRG_Protein_Result_rep2 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 2),]
HepaRG_Protein_Result_rep3 <-  HepaRG_Protein_Result[which(HepaRG_Protein_Result$biol_rep == 3),]

RPTEC_Protein_Result_rep1 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 1),]
RPTEC_Protein_Result_rep2 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 2),]
RPTEC_Protein_Result_rep3 <-  RPTEC_Protein_Result[which(RPTEC_Protein_Result$biol_rep == 3),]

## change rownames to apply merge
### function: enter dataframe name
### creates empty vector from rownames
### creates names for column 1-3
### apply make.names to be added to df
### add to df

change_rownames <- function(mydataframe){
  rowname_vector <- numeric(0)
  for(i in 1:nrow(mydataframe)){
    rowname <- paste(mydataframe[i,5], mydataframe[i,7], mydataframe[i,2], sep = "_")
    rowname_vector <- c(rowname_vector, rowname)
  }
  mydataframe <- cbind(rowname_vector, mydataframe)
  make.names(rowname_vector, unique = T)
  .rowNamesDF(x = mydataframe, make.names = T) <- rowname_vector
    return(mydataframe)
}

HepaRG_Protein_Result_rep1 <- change_rownames(HepaRG_Protein_Result_rep1)
HepaRG_Protein_Result_rep2 <- change_rownames(HepaRG_Protein_Result_rep2)
HepaRG_Protein_Result_rep3 <- change_rownames(HepaRG_Protein_Result_rep3)

RPTEC_Protein_Result_rep1 <- change_rownames(RPTEC_Protein_Result_rep1)
RPTEC_Protein_Result_rep2 <- change_rownames(RPTEC_Protein_Result_rep2)
RPTEC_Protein_Result_rep3 <- change_rownames(RPTEC_Protein_Result_rep3)

# calulate T/C per rep
##create function with input dataframe of rep
### attach df, create empty vectors
### go through whole df and look for substance_number which is not 0
### if so, go through whole df and look if a row with substance_number 0 matches regarding timepoint and concentration number
### if so, calculate T/C (100%), add to df
### df is list, has to be converted to list and colums converted to numeric
### return new df

calculate_T_C <- function(replicate){
  attach(replicate)
  newtable <- numeric(0)
  newthing <- numeric(0)
  for(i in 1:nrow(replicate)){
    if(substance_number[i] != 0){
      for(k in 1:nrow(replicate)){
        if(time_point[k] == time_point[i] & concentration_number[k] == concentration_number[i] & substance_number[k] == 0){
          newthing <- replicate[i, c(9:16)]/replicate[k,c(9:16)]*100
          rowofnewthing <- c(replicate[i,1:8], newthing)
          newtable <- rbind(newtable, rowofnewthing)
        }
      }
    }
  }
  detach(replicate)
  for(m in c(1,3,6)){
    newtable[,m] <- as.character(newtable[,m])
  }
  newname <- as.data.frame(newtable)
  rownames(newname) <- newname$rowname_vector
    for(l in c(2,4,5,7,8:16)){
    newname[,l] <- as.numeric(newname[,l])
  }
  return(newname)
}

HepaRG_Protein_Result_rep1_T_C <- calculate_T_C(HepaRG_Protein_Result_rep1)
HepaRG_Protein_Result_rep2_T_C <- calculate_T_C(HepaRG_Protein_Result_rep2)
HepaRG_Protein_Result_rep3_T_C <- calculate_T_C(HepaRG_Protein_Result_rep3)

RPTEC_Protein_Result_rep1_T_C <- calculate_T_C(RPTEC_Protein_Result_rep1)
RPTEC_Protein_Result_rep2_T_C <- calculate_T_C(RPTEC_Protein_Result_rep2)
RPTEC_Protein_Result_rep3_T_C <- calculate_T_C(RPTEC_Protein_Result_rep3)

# merge 3 dataframes
## first merge 2, merge third to that
## change rownames back, as they get set as new colomns
## get 3 rows with indications of the rows (as originally) from the new row.names
## create one dataframe and sort by substance name and concentration
## set rownames again
### input: dataframes, suffixes as "_suffix"

merge_3_data.frames <- function(dataframe1, suffix1, dataframe2, suffix2, dataframe3, suffix3){
  firstmerge <- merge(dataframe1[,c(9:16)],dataframe2[,c(9:16)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix1, suffix2))
  rownames(firstmerge) <- firstmerge$Row.names
  names(dataframe3)[9:16] <- paste(names(dataframe3)[9:16], "rep3", sep = "_")
  secondmerge <- merge(firstmerge[,c(2:ncol(firstmerge))], dataframe3[,c(9:ncol(dataframe3))], by = "row.names", all.x = T, all.y = T, suffixes = c("", suffix3))
  rownames(secondmerge) <- secondmerge$Row.names
  Beschriftung <- data.frame(t(data.frame(strsplit(as.character(secondmerge[c(1:length(secondmerge$Row.names)),1]), split = "_"))))
  dataframename <- cbind(Beschriftung, secondmerge[,c(2:ncol(secondmerge))])
  row.names(dataframename) <- secondmerge$Row.names
  names(dataframename)[1] <- "substance"
  names(dataframename)[2] <- "concentration_µM"
  names(dataframename)[3] <- "timepoint"
  return(dataframename)
}

HepaRG_Protein_Result_merge <- merge_3_data.frames(HepaRG_Protein_Result_rep1_T_C, "_rep1", HepaRG_Protein_Result_rep2_T_C, "_rep2", HepaRG_Protein_Result_rep3_T_C, "_rep3")

RPTEC_Protein_Result_merge <- merge_3_data.frames(RPTEC_Protein_Result_rep1_T_C, "_rep1", RPTEC_Protein_Result_rep2_T_C, "_rep2", RPTEC_Protein_Result_rep3_T_C, "_rep3")

#calculate mean
## create vector for new colnames
## create df with row information
### calculate rowmean for each analyte and attach it to df
### attach colnames
### same process for sd
## apply to merge df 
#statistics
## same process for stat
### extract p-value into a separate vector ---> add to df

install.packages("boot")
library(boot)

calculate_mean_sd <- function(dataframe_for_mean_sd){
  coloumnames <- gsub("_rep[0-9]", "", names(dataframe_for_mean_sd))
  coloumnames <- gsub("\\...*", "", coloumnames)
  mean_SD_df <- data.frame(dataframe_for_mean_sd[,c(1:3)])
  for (i in 4:11){
    analyte_mean <- rowMeans(dataframe_for_mean_sd[, c(seq(from = i, by = 8, length.out = 3))], na.rm = T)
    mean_SD_df <- cbind(mean_SD_df, analyte_mean)
    names(mean_SD_df)[i] <- paste(coloumnames[i], "mean", sep = "_")
  }
  for (i in 4:11){
    analyte_sd <- apply(dataframe_for_mean_sd[, c(seq(from = i, by = 8, length.out = 3))], MARGIN = 1, FUN = sd, na.rm = T)
    mean_SD_df <- cbind(mean_SD_df, analyte_sd)
    names(mean_SD_df)[i+8] <- paste(coloumnames[i], "sd", sep = "_")
  }
  statistic_function <- function(data, indices){
    sample_data <- data[indices]
    return(mean(sample_data))
  }
  for (i in 4:11) {
    t_test_substance <- numeric(0)
    p_value_as_asterisks <- numeric(0)
    bootstrap_as_asterisks <- numeric(0)
    for(k in 1:nrow(dataframe_for_mean_sd)){
      t_test <- t.test(as.numeric(dataframe_for_mean_sd[k, c(seq(from = i, by = 8, length.out = 3))]), mu = 100)
      t_test_substance <- c(t_test_substance, t_test$p.value)
      if(t_test$p.value < 0.05){
        as_asterisks <- "*"
      }else{
        as_asterisks <- ""
      }
      p_value_as_asterisks <- c(p_value_as_asterisks, as_asterisks)
      boot_result <- boot(na.omit(as.numeric(dataframe_for_mean_sd[k, c(seq(from = i, by = 8, length.out = 3))])), statistic_function, R = 1000)
      boostrap_means <- boot_result$t
      boot.ci(boot_result, type = "basic")
      confidence_interval <- boot.ci(boot_result, type = "basic")$basic
      if(confidence_interval[1,4] < 100 & 100 < confidence_interval[1,5]){
        bootstrap <- ""
      }else{
        bootstrap <- "*"
      }
      bootstrap_as_asterisks <- c(bootstrap_as_asterisks, bootstrap)
    }
    mean_SD_df <- cbind(mean_SD_df, t_test_substance, p_value_as_asterisks, bootstrap_as_asterisks)
    names(mean_SD_df) <- gsub("t_test_substance", paste(names(mean_SD_df[i]), "p-value", sep = "_"), names(mean_SD_df))
    names(mean_SD_df) <- gsub("p_value_as_asterisks", paste(names(mean_SD_df[i]), "asterisks", sep = "_"), names(mean_SD_df))
    names(mean_SD_df) <- gsub("bootstrap_as_asterisks", paste(names(mean_SD_df[i]), "boot", sep = "_"), names(mean_SD_df))
  }
  return(mean_SD_df)
}

HepaRG_Protein_Result_mean_SD <- calculate_mean_sd(HepaRG_Protein_Result_merge)

RPTEC_Protein_Result_mean_SD <- calculate_mean_sd(RPTEC_Protein_Result_merge)

HepaRG_Protein_Result_mean_SD <- cbind(row.names(HepaRG_Protein_Result_mean_SD), HepaRG_Protein_Result_mean_SD)
names(HepaRG_Protein_Result_mean_SD)[1] <- "longname"

RPTEC_Protein_Result_mean_SD <- cbind(row.names(RPTEC_Protein_Result_mean_SD), RPTEC_Protein_Result_mean_SD)
names(RPTEC_Protein_Result_mean_SD)[1] <- "longname"

library(ggplot2)
library(reshape)
library(tidyverse)

# heatmap with significance
## transform p-value to asterisks
## create melt with asterisks
## transform mean_values to factors
### can be included in heatmap --> creates discrete heatmap, one color for a certain range

data_mean_for_heatmap <- function(dataframe){
  datamelt_mean <- melt(dataframe[, c(1, 5:12)])
  datamelt_mean$variable <- gsub("_mean", "", datamelt_mean$variable)
  datamelt_mean$value_group <- cut(datamelt_mean$value, 
                                   breaks = c(0, 50, 75, 150, 200, Inf), 
                                   labels = c("< 50", "50-75", "75-150", "150-200", "> 200"),
                                   right = FALSE)
  datamelt_mean$longname <- gsub("_", " ", datamelt_mean$longname)
  datamelt_mean$longname <- gsub("\\s(36h|72h)", "µM \\1", datamelt_mean$longname, perl = TRUE)
  return(datamelt_mean)
}

data_asterisks_for_heatmap <- function(dataframe){
  datamelt_stat <- melt(dataframe[, c(1, grep("_boot", names(dataframe)))], id.vars = "longname")
  datamelt_stat$variable <- gsub("_mean_boot", "", datamelt_stat$variable)
  datamelt_stat$longname <- gsub("_", " ", datamelt_stat$longname)
  datamelt_stat$longname <- gsub("\\s(36h|72h)", "µM \\1", datamelt_stat$longname, perl = TRUE)
  return(datamelt_stat)
}

### black box around substances
lines_heatmap <- function(dataframe){
  boxy <- data.frame(xmin = 0.5,
                   xmax = 8.5,
                   ymin = seq(0.5, nrow(dataframe) - 0.5, 4),
                   ymax = seq(4.5, nrow(dataframe) + 0.5, 4))
  return(boxy)
}

datamelt_mean_HepaRG <- data_mean_for_heatmap(HepaRG_Protein_Result_mean_SD)
datamelt_stat_HepaRG <- data_asterisks_for_heatmap(HepaRG_Protein_Result_mean_SD)
boxy_HepaRG <- lines_heatmap(HepaRG_Protein_Result_mean_SD)

datamelt_mean_RPTEC <- data_mean_for_heatmap(RPTEC_Protein_Result_mean_SD)
datamelt_stat_RPTEC <- data_asterisks_for_heatmap(RPTEC_Protein_Result_mean_SD)
boxy_RPTEC <- lines_heatmap(RPTEC_Protein_Result_mean_SD)

### add space between number and unit

add_spaces <- function(text) {
  text <- gsub("(\\d)µM", "\\1 µM", text)
  text <- gsub("(\\d)h", "\\1 h", text)
  return(text)
}

datamelt_mean_HepaRG$longname <- add_spaces(datamelt_mean_HepaRG$longname)
datamelt_stat_HepaRG$longname <- add_spaces(datamelt_stat_HepaRG$longname)

datamelt_mean_RPTEC$longname <- add_spaces(datamelt_mean_RPTEC$longname)
datamelt_stat_RPTEC$longname <-add_spaces(datamelt_stat_RPTEC$longname)

### change reihenfolge der y-achse

datamelt_mean_HepaRG$longname <- as.factor(datamelt_mean_HepaRG$longname)
datamelt_mean_RPTEC$longname <- as.factor(datamelt_mean_RPTEC$longname)
datamelt_stat_HepaRG$longname <- as.factor(datamelt_stat_HepaRG$longname)
datamelt_stat_RPTEC$longname <- as.factor(datamelt_stat_RPTEC$longname)
orderHepaRG <- c(2, 1, 4, 3, 24, 23, 22, 21, 12, 11, 10, 9, 8, 7, 6, 5, 20, 19, 18, 17, 14, 13, 16, 15) 
orderRPTEC <- c(2, 1, 4, 3, 24, 23, 22, 21, 12, 11, 10, 9, 8, 7, 6, 5, 20, 19, 18, 17, 16, 15, 14, 13)
datamelt_stat_HepaRG$longname <- factor(datamelt_stat_HepaRG$longname, levels = levels(datamelt_mean_HepaRG$longname)[orderHepaRG])
datamelt_mean_HepaRG$longname <- factor(datamelt_mean_HepaRG$longname,  levels = levels(datamelt_mean_HepaRG$longname)[orderHepaRG])
datamelt_stat_RPTEC$longname <- factor(datamelt_stat_RPTEC$longname, levels = levels(datamelt_mean_RPTEC$longname)[orderRPTEC])
datamelt_mean_RPTEC$longname <- factor(datamelt_mean_RPTEC$longname,  levels = levels(datamelt_mean_RPTEC$longname)[orderRPTEC])

# heatmap
## data is displayed with numbers (geom_text)
## discrete colors (see cut function above)
## geom_rect for lines
## geom_text(2) for asterisks in upper right corner

library(grDevices)

###RPTEC

Heatmap_RPTEC <- ggplot() +
  geom_raster(data = datamelt_mean_RPTEC, aes(x = variable, y = longname, fill = value_group), show.legend = T) +
  scale_fill_manual(name = "fold change (%)", 
                    values = c("< 50" = "dodgerblue", "50-75" = "lightskyblue", "75-150" = "white", "150-200" = "lightpink", "> 200" = "red2"),
                    drop = F,
                    ) +
  geom_rect(data = boxy_RPTEC, 
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax), 
            colour = "black", fill = NA, size = 1)+
  geom_text(data = datamelt_mean_RPTEC, aes(x = variable, y = longname, label = round(value)), size = 3) +
  geom_text(data = datamelt_stat_RPTEC, aes(x = variable, y = longname, label = value), 
            nudge_y = 0.05, nudge_x = 0.25, size = 5)+
  theme(axis.title = element_blank(), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10))

png(filename = "RPTEC_Luminex.png",
    width = 210,
    height = 180,
    units = "mm",
    bg = "white",
    res = 900)

Heatmap_RPTEC

dev.off()

###HepaRG

Heatmap_HepaRG <- ggplot() +
  geom_raster(data = datamelt_mean_HepaRG, aes(x = variable, y = longname, fill = value_group)) +
  scale_fill_manual(name = "fold change (%)", 
                    values = c("dodgerblue", "lightskyblue1", "white", "lightpink", "red2"),
                    drop = FALSE) +
  geom_rect(data = boxy_HepaRG, 
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax), 
            colour = "black", fill = NA, size = 1)+
  geom_text(data = datamelt_mean_HepaRG, aes(x = variable, y = longname, label = round(value)), size = 3) +
  geom_text(data = datamelt_stat_HepaRG, aes(x = variable, y = longname, label = value), 
            nudge_y = 0.05, nudge_x = 0.25, size = 5)+
  theme(axis.title = element_blank(), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), legend.position = "right")


png(filename = "HepaRG_Luminex.png",
    width = 210,
    height = 180,
    units = "mm",
    bg = "white",
    res = 900)

Heatmap_HepaRG

dev.off()

#spupplementary material
##boxplot
###HepaRG

HepaRG_Protein_datamelt_boxplot <- melt(cbind(longnames = row.names(HepaRG_Protein_Result_merge), HepaRG_Protein_Result_merge))

HepaRG_Protein_datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", HepaRG_Protein_datamelt_boxplot$variable)

Sup_Fig_HepaRG <- ggplot(HepaRG_Protein_datamelt_boxplot, aes(longnames, value))+
  geom_point(size = 1)+
  geom_boxplot(width = 0.5, size = 0.3)+
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.2,
    size = 0.3
  )+
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1
  )+
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 0.5)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank())+
  ylab("relative protein expression T/C (%)")+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  ggtitle("HepaRG")

ggsave("Sup_Fig_HepaRG.png",width = 9.21, height = 6.85, dpi = 1200)
 
###RPTEC

RPTEC_Protein_datamelt_boxplot <- melt(cbind(longnames = row.names(RPTEC_Protein_Result_merge), RPTEC_Protein_Result_merge))

RPTEC_Protein_datamelt_boxplot$variable <- gsub("_rep[0-9]+", "", RPTEC_Protein_datamelt_boxplot$variable)

ggplot(RPTEC_Protein_datamelt_boxplot, aes(longnames, value))+
  geom_point(size = 1)+
  geom_boxplot(width = 0.5, size = 0.3)+
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.2,
    size = 0.3
  )+
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1
  )+
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 0.5)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank())+
  ylab("relative protein expression T/C (%)")+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  ggtitle("RPTEC")

ggsave("Sup_Fig_RPTEC_Luminex.png",width = 9.21, height = 6.85, dpi = 1200)



