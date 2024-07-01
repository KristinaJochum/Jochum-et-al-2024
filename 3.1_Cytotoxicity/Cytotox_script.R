#set: import data and names, Barchart name
#check: SC position

#import all data

## data files as csv

### get library
library(readr)

### generate list of all files with .csv (=datafiles)
temp <- list.files(pattern = "*.csv")

### read all files in list from above, generate name by removing .csv, assign file with corresponding filename, repeat till end of list
for (i in 1:length(temp)) {
  file <- read_delim(temp[i], delim = ";", escape_double = FALSE, trim_ws = TRUE)
  filename <- make.names(gsub(".csv", "", temp[i]))
  assign(filename, file)
}

## name files as excel to get µM

### get library
library(readxl)

### generate a list with all .xlsx files(=namefiles)
temp_names <- list.files(pattern="*.xlsx")

### read all files in list from above, generate name by removing .xlsx, assign file with corresponding filename, repeat till end of list
for (i in 1:length(temp_names)) {
  file <- read_excel(temp_names[i], na = "Null")
  filename <- make.names(gsub(".xlsx", "", temp_names[i]))
  assign(filename, file)
}

### remove variables that are not needed
rm(file,filename,i,temp,temp_names)


# dataframe for 1 Experiment
## create function
### subtract blank from each well (blank is above/beneave triplicate)
### if blank is in upper half of plate, select "oben" as blankposition, otherwise select "unten"
### labeling of data is generated from names file (mynames), split characterstring to be able to sort by concentration and name
### put everxthing in one barchart
### save barchart in corresponding directory
### print barchart and check SC (should be 100)

generate_dataframe <- function(mydata, mynames, blank, blankposition){
  minusBlank_oben <-  colMeans(mydata[c(2,3,4), c(seq(2,11,1))] - mydata[c(1,1,1), c(seq(2,11,1))])
  minusBlank_unten <-  colMeans(mydata[c(5,6,7), c(seq(2,11,1))] - mydata[c(8,8,8), c(seq(2,11,1))])
  if (blankposition == "oben"){
    T_C_oben <- minusBlank_oben*100/minusBlank_oben[blank]
    T_C_unten <- minusBlank_unten*100/minusBlank_oben[blank]
    T_C <- c(T_C_oben, T_C_unten)
    SD_unten <- apply(mydata[c(5,6,7), c(seq(2,11,1))] - mydata[c(8,8,8), c(seq(2,11,1))], MARGIN = 2,FUN = "sd")*100/minusBlank_oben[blank]
    SD_oben <- apply(mydata[c(2,3,4), c(seq(2,11,1))] - mydata[c(1,1,1), c(seq(2,11,1))],MARGIN = 2,FUN = "sd")*100/minusBlank_oben[blank]
    SD <- c(SD_oben, SD_unten)
  } else {
    T_C_oben <- minusBlank_oben*100/minusBlank_unten[blank]
    T_C_unten <- minusBlank_unten*100/minusBlank_unten[blank]
    T_C <- c(T_C_oben, T_C_unten)
    SD_unten <- apply(mydata[c(5,6,7), c(seq(2,11,1))] - mydata[c(8,8,8), c(seq(2,11,1))], MARGIN = 2,FUN = "sd")*100/minusBlank_unten[blank]
    SD_oben <- apply(mydata[c(2,3,4), c(seq(2,11,1))] - mydata[c(1,1,1), c(seq(2,11,1))],MARGIN = 2,FUN = "sd")*100/minusBlank_unten[blank]
    SD <- c(SD_oben, SD_unten)
  }
  Beschriftung <- as.data.frame(t(as.data.frame(strsplit(unlist(c(mynames[1, c(seq(2,11,1))], mynames[8, c(seq(2,11,1))])), split = " "))))
  Beschriftung[,2] <- as.numeric(Beschriftung[,2])
  DF <- data.frame(Beschriftung, T_C, SD)
  DF <- DF[order(DF$V1, DF$V2),]
  return(DF)
}

## include in generate_Barchart
### mydata: datafile 
### mynames: namefile
### blank: blankposition in "" indicated by colname (equivalent to position on wellplate) not actual position in the dataframe

RPTEC_rep3_1_3_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep3_1_3_72h_WST1, RPTEC_rep3_1_3_72h_names, "8", "oben")
RPTEC_rep3_1_3_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep3_1_3_72h_NRU, RPTEC_rep3_1_3_72h_names, "8", "oben")

RPTEC_rep3_4_6_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep3_4_6_72h_WST1, RPTEC_rep3_4_6_72h_names, "8", "oben")
RPTEC_rep3_4_6_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep3_4_6_72h_NRU, RPTEC_rep3_4_6_72h_names, "8", "oben")

RPTEC_rep4_1_3_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep4_1_3_72h_WST1, RPTEC_rep4_1_3_72h_names, "8", "oben")

RPTEC_rep4_4_6_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep4_4_6_72h_WST1, RPTEC_rep4_4_6_72h_names, "8", "oben")

RPTEC_rep5_1_3_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep5_1_3_72h_WST1, RPTEC_rep5_1_3_72h_names, "8", "oben")
RPTEC_rep5_1_3_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep5_1_3_72h_NRU, RPTEC_rep5_1_3_72h_names, "8", "oben")

RPTEC_rep5_4_6_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep5_4_6_72h_WST1, RPTEC_rep5_4_6_72h_names, "8", "oben")
RPTEC_rep5_4_6_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep5_4_6_72h_NRU, RPTEC_rep5_4_6_72h_names, "8", "oben")

RPTEC_rep6_1_3_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep6_1_3_72h_WST1, RPTEC_rep6_1_3_72h_names, "8", "oben")
RPTEC_rep6_1_3_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep6_1_3_72h_NRU, RPTEC_rep6_1_3_72h_names, "8", "oben")

RPTEC_rep6_4_6_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep6_4_6_72h_WST1, RPTEC_rep6_4_6_72h_names, "8", "oben")
RPTEC_rep6_4_6_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep6_4_6_72h_NRU, RPTEC_rep6_4_6_72h_names, "8", "oben")

RPTEC_rep7_1_2_6_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep7_1_2_6_72h_WST1, RPTEC_rep7_1_2_6_72h_names, "8", "oben")
RPTEC_rep7_1_2_6_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep7_1_2_6_72h_NRU, RPTEC_rep7_1_2_6_72h_names, "8", "oben")

RPTEC_rep8_1_3_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep8_1_3_72h_WST1, RPTEC_rep8_1_3_72h_names, "8", "oben")
RPTEC_rep8_1_3_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep8_1_3_72h_NRU, RPTEC_rep8_1_3_72h_names, "8", "oben")

RPTEC_rep8_4_6_72h_WST1_dataframe <- generate_dataframe(RPTEC_rep8_4_6_72h_WST1, RPTEC_rep8_4_6_72h_names, "8", "oben")
RPTEC_rep8_4_6_72h_NRU_dataframe <- generate_dataframe(RPTEC_rep8_4_6_72h_NRU, RPTEC_rep8_4_6_72h_names, "8", "oben")

#create one dataframe per replicate_WST1
##rep3
RPTEC_rep3_72h_WST1_dataframe <- rbind(
  RPTEC_rep3_1_3_72h_WST1_dataframe, 
  RPTEC_rep3_4_6_72h_WST1_dataframe)

### sort
RPTEC_rep3_72h_WST1_dataframe <- RPTEC_rep3_72h_WST1_dataframe[order(RPTEC_rep3_72h_WST1_dataframe$V1, RPTEC_rep3_72h_WST1_dataframe$V2),]

## rep4
RPTEC_rep4_72h_WST1_dataframe <- rbind(
  RPTEC_rep4_1_3_72h_WST1_dataframe, 
  RPTEC_rep4_4_6_72h_WST1_dataframe)

### sort
RPTEC_rep4_72h_WST1_dataframe <- RPTEC_rep4_72h_WST1_dataframe[order(RPTEC_rep4_72h_WST1_dataframe$V1, RPTEC_rep4_72h_WST1_dataframe$V2),]

## rep5
RPTEC_rep5_72h_WST1_dataframe <- rbind(
  RPTEC_rep5_1_3_72h_WST1_dataframe, 
  RPTEC_rep5_4_6_72h_WST1_dataframe)

### sort
RPTEC_rep5_72h_WST1_dataframe <- RPTEC_rep5_72h_WST1_dataframe[order(RPTEC_rep5_72h_WST1_dataframe$V1, RPTEC_rep5_72h_WST1_dataframe$V2),]

## rep6
RPTEC_rep6_72h_WST1_dataframe <- rbind(
  RPTEC_rep6_1_3_72h_WST1_dataframe, 
  RPTEC_rep6_4_6_72h_WST1_dataframe)

### sort
RPTEC_rep6_72h_WST1_dataframe <- RPTEC_rep6_72h_WST1_dataframe[order(RPTEC_rep6_72h_WST1_dataframe$V1, RPTEC_rep6_72h_WST1_dataframe$V2),]

## rep7
RPTEC_rep7_72h_WST1_dataframe <- RPTEC_rep7_1_2_6_72h_WST1_dataframe

### sort
RPTEC_rep7_72h_WST1_dataframe <- RPTEC_rep7_72h_WST1_dataframe[order(RPTEC_rep7_72h_WST1_dataframe$V1, RPTEC_rep7_72h_WST1_dataframe$V2),]


## rep8
RPTEC_rep8_72h_WST1_dataframe <- rbind(
  RPTEC_rep8_1_3_72h_WST1_dataframe, 
  RPTEC_rep8_4_6_72h_WST1_dataframe)

### sort
RPTEC_rep8_72h_WST1_dataframe <- RPTEC_rep8_72h_WST1_dataframe[order(RPTEC_rep8_72h_WST1_dataframe$V1, RPTEC_rep8_72h_WST1_dataframe$V2),]

#create one dataframe per replicate_NRU
##rep3
RPTEC_rep3_72h_NRU_dataframe <- rbind(
  RPTEC_rep3_1_3_72h_NRU_dataframe, 
  RPTEC_rep3_4_6_72h_NRU_dataframe)

### sort
RPTEC_rep3_72h_NRU_dataframe <- RPTEC_rep3_72h_NRU_dataframe[order(RPTEC_rep3_72h_NRU_dataframe$V1, RPTEC_rep3_72h_NRU_dataframe$V2),]

## rep5
RPTEC_rep5_72h_NRU_dataframe <- rbind(
  RPTEC_rep5_1_3_72h_NRU_dataframe, 
  RPTEC_rep5_4_6_72h_NRU_dataframe)

### sort
RPTEC_rep5_72h_NRU_dataframe <- RPTEC_rep5_72h_NRU_dataframe[order(RPTEC_rep5_72h_NRU_dataframe$V1, RPTEC_rep5_72h_NRU_dataframe$V2),]

## rep6
RPTEC_rep6_72h_NRU_dataframe <- rbind(
  RPTEC_rep6_1_3_72h_NRU_dataframe, 
  RPTEC_rep6_4_6_72h_NRU_dataframe)

### sort
RPTEC_rep6_72h_NRU_dataframe <- RPTEC_rep6_72h_NRU_dataframe[order(RPTEC_rep6_72h_NRU_dataframe$V1, RPTEC_rep6_72h_NRU_dataframe$V2),]

## rep7
RPTEC_rep7_72h_NRU_dataframe <- RPTEC_rep7_1_2_6_72h_NRU_dataframe

### sort
RPTEC_rep7_72h_NRU_dataframe <- RPTEC_rep7_72h_NRU_dataframe[order(RPTEC_rep7_72h_NRU_dataframe$V1, RPTEC_rep7_72h_NRU_dataframe$V2),]


## rep8
RPTEC_rep8_72h_NRU_dataframe <- rbind(
  RPTEC_rep8_1_3_72h_NRU_dataframe, 
  RPTEC_rep8_4_6_72h_NRU_dataframe)

### sort
RPTEC_rep8_72h_NRU_dataframe <- RPTEC_rep8_72h_NRU_dataframe[order(RPTEC_rep8_72h_NRU_dataframe$V1, RPTEC_rep8_72h_NRU_dataframe$V2),]


# add all dataframes together
## change rownames to apply merge
### function: enter barchart name
### creates empty vector for rownames
### creates names for column 1-3
### apply make.names to be added to df
### add to df

change_rownames <- function(mydataframe){
  rowname_vector <- numeric(0)
  for(i in 1:length(mydataframe[,1])){
    rowname <- paste(mydataframe[i,1], mydataframe[i,2], mydataframe[i,3], sep = "_")
    rowname_vector <- c(rowname_vector, rowname)
    rowname_vector <- gsub("2-", "X2.", rowname_vector)
  }
  make.names(rowname_vector, unique = T)
  .rowNamesDF(x = mydataframe, make.names = T) <- rowname_vector
  return(mydataframe)
}

RPTEC_rep3_72h_WST1_dataframe <- change_rownames(RPTEC_rep3_72h_WST1_dataframe)

RPTEC_rep4_72h_WST1_dataframe <- change_rownames(RPTEC_rep4_72h_WST1_dataframe)

RPTEC_rep5_72h_WST1_dataframe <- change_rownames(RPTEC_rep5_72h_WST1_dataframe)

RPTEC_rep6_72h_WST1_dataframe <- change_rownames(RPTEC_rep6_72h_WST1_dataframe)

RPTEC_rep7_72h_WST1_dataframe <- change_rownames(RPTEC_rep7_72h_WST1_dataframe)

RPTEC_rep8_72h_WST1_dataframe <- change_rownames(RPTEC_rep8_72h_WST1_dataframe)


RPTEC_rep3_72h_NRU_dataframe <- change_rownames(RPTEC_rep3_72h_NRU_dataframe)

RPTEC_rep5_72h_NRU_dataframe <- change_rownames(RPTEC_rep5_72h_NRU_dataframe)

RPTEC_rep6_72h_NRU_dataframe <- change_rownames(RPTEC_rep6_72h_NRU_dataframe)

RPTEC_rep7_72h_NRU_dataframe <- change_rownames(RPTEC_rep7_72h_NRU_dataframe)

RPTEC_rep8_72h_NRU_dataframe <- change_rownames(RPTEC_rep8_72h_NRU_dataframe)

# merge 5 dataframes
## first merge 2, then merge these 2, merge fifth to that
## change rownames back, as they get set as new colomns
## get 3 rows with indications of the rows (as originally) from the new row.names
## create one dataframe and sort by substance name and concentration
## set rownames again
### input: dataframes, suffixes as "_suffix"

merge_5_data.frames <- function(dataframe1, suffix1, dataframe2, suffix2, dataframe3, suffix3, dataframe4, suffix4, dataframe5, suffix5){
  firstmerge <- merge(dataframe1[,c(4,5)],dataframe2[,c(4,5)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix1, suffix2))
  rownames(firstmerge) <- firstmerge$Row.names
  secondmerge <- merge(dataframe3[,c(4,5)],dataframe4[,c(4,5)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix3, suffix4))
  rownames(secondmerge) <- secondmerge$Row.names
  thirdmerge <- merge(firstmerge[,c(2,3,4,5)],secondmerge[,c(2:5)], by = "row.names", all.x = T, all.y = T)
  rownames(thirdmerge) <- thirdmerge$Row.names
  forthmerge <- merge(thirdmerge[,c(2:9)], dataframe5[,c(4,5)], by = "row.names", all.x = T, all.y = T, suffixes = c("", suffix5))
  rownames(forthmerge) <- forthmerge$Row.names
  Beschriftung <- data.frame(t(data.frame(strsplit(as.character(forthmerge[c(1:length(forthmerge$Row.names)),1]), split = "_"))))
  Beschriftung[,2] <- as.numeric(Beschriftung[,2])
  dataframename <- cbind(Beschriftung, forthmerge[,c(2:length(forthmerge))])
  row.names(dataframename) <- forthmerge[,1]
  dataframename <- dataframename[order(dataframename[,1], dataframename[,2]),]
  names(dataframename)[1] <- "substance"
  names(dataframename)[2] <- "concentration"
  names(dataframename)[3] <- "unit"
  View(dataframename)
  return(dataframename)
}

merge_6_data.frames <- function(dataframe1, suffix1, dataframe2, suffix2, dataframe3, suffix3, dataframe4, suffix4, dataframe5, suffix5, dataframe6, suffix6){
  print(paste0(suffix1, suffix2, suffix3, suffix4, suffix5, suffix6))
  firstmerge <- merge(dataframe1[,c(4,5)], dataframe2[,c(4,5)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix1, suffix2))
  rownames(firstmerge) <- firstmerge$Row.names
  secondmerge <- merge(dataframe3[,c(4,5)], dataframe4[,c(4,5)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix3, suffix4))
  rownames(secondmerge) <- secondmerge$Row.names
  thirdmerge <- merge(firstmerge[,c(2,3,4,5)],secondmerge[,c(2:5)], by = "row.names", all.x = T, all.y = T)
  rownames(thirdmerge) <- thirdmerge$Row.names
  forthmerge <- merge(dataframe5[,c(4,5)], dataframe6[,c(4,5)], by = "row.names", all.x = T, all.y = T, suffixes = c(suffix5, suffix6))
  rownames(forthmerge) <- forthmerge$Row.names
  fifthmerge <- merge(thirdmerge[,c(2:9)], forthmerge[,c(2:5)], by = "row.names", all.x = T, all.y = T)
  rownames(fifthmerge) <- fifthmerge$Row.names
  Beschriftung <- data.frame(t(data.frame(strsplit(as.character(fifthmerge[c(1:length(fifthmerge$Row.names)),1]), split = "_"))))
  Beschriftung[,2] <- as.numeric(Beschriftung[,2])
  dataframename <- cbind(Beschriftung, fifthmerge[,c(2:length(fifthmerge))])
  row.names(dataframename) <- fifthmerge[,1]
  dataframename <- dataframename[order(dataframename[,1], dataframename[,2]),]
  names(dataframename)[1] <- "substance"
  names(dataframename)[2] <- "concentration"
  names(dataframename)[3] <- "unit"
  return(dataframename)
}

RPTEC_allrep_72h_WST1_dataframe <- merge_6_data.frames(RPTEC_rep3_72h_WST1_dataframe, "_rep3", RPTEC_rep4_72h_WST1_dataframe, "_rep4", RPTEC_rep5_72h_WST1_dataframe, "_rep5", RPTEC_rep6_72h_WST1_dataframe, "_rep6", RPTEC_rep7_72h_WST1_dataframe, "_rep7", RPTEC_rep8_72h_WST1_dataframe, "_rep8")
RPTEC_allrep_72h_NRU_dataframe <- merge_5_data.frames(RPTEC_rep3_72h_NRU_dataframe, "_rep3", RPTEC_rep5_72h_NRU_dataframe, "_rep5", RPTEC_rep6_72h_NRU_dataframe, "_rep6", RPTEC_rep7_72h_NRU_dataframe, "_rep7", RPTEC_rep8_72h_NRU_dataframe, "_rep8")


library(ggplot2)

#create dataframe with means and SD, substances only
##substances only (rownumber)

Cytotox_mean_SD <- function(dataframe){
  dataframe$substance <- gsub("X2.P", "2-P", dataframe$substance)
  needed_rows <- numeric()
  pesticides <- c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "X2.Phenylphenol")
  n_above_3 <- which(rowSums(!is.na(dataframe[, c(seq(4, ncol(dataframe), 2))])) >= 3)
  print(names(n_above_3))
  for(i in pesticides){
    row_addition <- n_above_3[grep(pattern = i, names(n_above_3))]
    needed_rows <- c(needed_rows, row_addition)
  }
  substance_means <- rowMeans(dataframe[needed_rows,c(seq(4, ncol(dataframe), 2))], na.rm = T)
  substance_SD <- apply(X = dataframe[needed_rows,c(seq(4, ncol(dataframe), 2))], MARGIN = 1, FUN = sd, na.rm = T)
  means_SD_allrep_barchart <- data.frame(
    substance = factor(dataframe$substance[needed_rows], c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "2-Phenylphenol")), 
    concentration = dataframe$concentration[needed_rows], 
    unit = dataframe$unit[needed_rows], 
    substance_means,
    substance_SD,
    row.names = row.names(dataframe[needed_rows,])
  )
  return(means_SD_allrep_barchart)
}

RPTEC_means_SD_allrep_WST1_barchart <- Cytotox_mean_SD(RPTEC_allrep_72h_WST1_dataframe)
RPTEC_means_SD_allrep_NRU_barchart <- Cytotox_mean_SD(RPTEC_allrep_72h_NRU_dataframe)

#create basic bar plot
##define data, concentration as factor
df_base_WST1 <- ggplot(
  data = RPTEC_means_SD_allrep_WST1_barchart,
  aes(x = factor(concentration), y = ifelse(RPTEC_means_SD_allrep_WST1_barchart$substance_means < 0, 0, RPTEC_means_SD_allrep_WST1_barchart$substance_means))
)

##specify tings
### errrorbar width and everything, horizontal line, bar width, facetwrap, theme
RPTEC_Cytotox_WST1 <- df_base_WST1 + 
  geom_errorbar(aes(
    ymin = ifelse(RPTEC_means_SD_allrep_WST1_barchart$substance_means -
      RPTEC_means_SD_allrep_WST1_barchart$substance_SD < 0, 0, RPTEC_means_SD_allrep_WST1_barchart$substance_means -
        RPTEC_means_SD_allrep_WST1_barchart$substance_SD) ,
    ymax = ifelse(RPTEC_means_SD_allrep_WST1_barchart$substance_means +
      RPTEC_means_SD_allrep_WST1_barchart$substance_SD < 0, 0, RPTEC_means_SD_allrep_WST1_barchart$substance_means +
        RPTEC_means_SD_allrep_WST1_barchart$substance_SD)
  ), width = 0.2, linewidth = 0.3)+
  geom_hline(aes(yintercept = 100), linetype = 2, size = 0.5, linewidth = 0.3)+
  geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.3)+
  facet_wrap(~substance, scales = "free_x", ncol = 6)+
  labs(title = "a", x = "concentration (µM)", y = "cell viability T/C (%)")+
  theme_bw()+
  theme(axis.text = element_text(size = 6),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        plot.title.position = "plot",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 235), breaks = c(seq(50, 200, 50)))

RPTEC_Cytotox_WST1

#create basic bar plot
##define data, concentration as factor
df_base_NRU <- ggplot(
  data = RPTEC_means_SD_allrep_NRU_barchart,
  aes(x = factor(concentration), y = ifelse(RPTEC_means_SD_allrep_NRU_barchart$substance_means < 0, 0, RPTEC_means_SD_allrep_NRU_barchart$substance_means))
)

##specify tings
### errrorbar width and everything, horizontal line, bar width, facetwrap, theme
RPTEC_Cytotox_NRU <- df_base_NRU + 
  geom_errorbar(aes(
    ymin = ifelse(RPTEC_means_SD_allrep_NRU_barchart$substance_means -
                    RPTEC_means_SD_allrep_NRU_barchart$substance_SD < 0, 0, RPTEC_means_SD_allrep_NRU_barchart$substance_means -
                    RPTEC_means_SD_allrep_NRU_barchart$substance_SD) ,
    ymax = ifelse(RPTEC_means_SD_allrep_NRU_barchart$substance_means +
                    RPTEC_means_SD_allrep_NRU_barchart$substance_SD < 0, 0, RPTEC_means_SD_allrep_NRU_barchart$substance_means +
                    RPTEC_means_SD_allrep_NRU_barchart$substance_SD)
  ), width = 0.2, linewidth = 0.3)+
  geom_hline(aes(yintercept = 100), linetype = 2, linewidth = 0.3)+
  geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.3)+
  facet_wrap(~substance, scales = "free_x", ncol = 6)+
  labs(title = "b", x = "concentration (µM)", y = "cell viability T/C (%)")+
  theme_bw()+
  theme(axis.text = element_text(size = 6),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        plot.title.position = "plot",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 10))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 235), breaks = c(seq(50, 200, 50)))

RPTEC_Cytotox_NRU

library(gridExtra)

RPTEC_Cytotox <- grid.arrange(RPTEC_Cytotox_WST1, RPTEC_Cytotox_NRU, ncol = 1)

ggsave(
  filename = paste0("Cytotox_", ".png"),
  plot = RPTEC_Cytotox,
  device = "png",
  width = 220,
  height = 130,
  units = "mm",
  dpi = 900)






