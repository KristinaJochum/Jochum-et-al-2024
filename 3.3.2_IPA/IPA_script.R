#IPA Graph for Paper
#Y-Axis: Diseases (top 10 overall), X-Axis: p-value, Bars: HepaRG, RPTEC, combined

#1.upload all results from HepaRG, RPTEC, combi
#2.combine all of theme in one df by effect
#3.determine the effects coming up the most based on p-value as well
#4.make bar graph as described above

## get library
library(readr)
library(tidyverse)

## data files as txt
### generate list of all files .txt (=datafiles)
### read all files from list, 
### skip first 2 rows, only cat and p-value cols
### remove doubled categories and only keep the one with the lowest p-value
#### Group by category and select the row with the minimum number,  ungroup
### change colname to include cell line, generate name from temp as substance
### assign file with corresponding filename, and add to result_list

Import <- function(CellLine){
  setwd(paste0("IPA_", CellLine, "/"))
  temp <- list.files(pattern = "*.txt")
  result_list <- list()
  for (i in 1:length(temp)) {
    file <- read_delim(temp[i], delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 2, locale = locale(decimal_mark = ","))
    file <- file[,c(1,3)]
    result_df <- file %>%
      group_by(Categories) %>%
      summarise(`p-value` = min(`p-value`))
    result_df <- result_df %>% ungroup()
    colnames(result_df) <- gsub("p-value", paste(colnames(result_df[2]), CellLine, gsub(".txt", "", temp[i]), sep = "_"), colnames(result_df), )
    filename <- paste(gsub(".txt", "", temp[i]))
    result_list[[filename]] <- result_df
  }
  setwd(("../"))
  return(result_list)
}

KSJ_HepaRG_IPA <- Import("HepaRG")
KSJ_RPTEC_IPA <- Import("RPTEC")
KSJ_combined_IPA <- Import("combined")

### combine more than one df
#### generate list of all df

df_list <- list(as.data.frame(KSJ_HepaRG_IPA[["Cyproconazole"]]),
                as.data.frame(KSJ_HepaRG_IPA[["Fluxapyroxad"]]), 
                as.data.frame(KSJ_HepaRG_IPA[["Azoxystrobin"]]),
                as.data.frame(KSJ_HepaRG_IPA[["Chlorotoluron"]]),
                as.data.frame(KSJ_HepaRG_IPA[["Thiabendazole"]]),
                as.data.frame(KSJ_HepaRG_IPA[["2Phenylphenol"]]),
                as.data.frame(KSJ_RPTEC_IPA[["Cyproconazole"]]),
                as.data.frame(KSJ_RPTEC_IPA[["Fluxapyroxad"]]), 
                as.data.frame(KSJ_RPTEC_IPA[["Azoxystrobin"]]),
                as.data.frame(KSJ_RPTEC_IPA[["Chlorotoluron"]]),
                as.data.frame(KSJ_RPTEC_IPA[["Thiabendazole"]]),
                as.data.frame(KSJ_RPTEC_IPA[["2Phenylphenol"]]),
                as.data.frame(KSJ_combined_IPA[["Cyproconazole"]]),
                as.data.frame(KSJ_combined_IPA[["Fluxapyroxad"]]), 
                as.data.frame(KSJ_combined_IPA[["Azoxystrobin"]]),
                as.data.frame(KSJ_combined_IPA[["Chlorotoluron"]]),
                as.data.frame(KSJ_combined_IPA[["Thiabendazole"]]),
                as.data.frame(KSJ_combined_IPA[["2Phenylphenol"]])
                )

#### merged by Categories

allCat_pvalue <- Reduce(function(x, y) merge(x, y, by = "Categories", all=TRUE), df_list)

###remove non-liver/-kidney effects

Cardiac_effects <- c(
  grep("Heart", allCat_pvalue$Categories),
  grep("Cardiac", allCat_pvalue$Categories),
  grep("Pulmonary", allCat_pvalue$Categories)
)

Cat_pvalue <- allCat_pvalue[-Cardiac_effects, ]
row.names(Cat_pvalue) <- Cat_pvalue$Categories

Cat_pvalue$rowSums <- rowSums(!is.na(Cat_pvalue))

top10 <- tail(Cat_pvalue[order(rowSums(!is.na(Cat_pvalue))),], 10)
library(reshape)
top10_melt <- melt(top10[, c(1:19)])
Beschriftung <- data.frame(t(data.frame(strsplit(as.character(top10_melt[,2]), split = "_"))))
top10_melt <- cbind(Beschriftung[, c(2,3)], top10_melt)

##graph
library(ggplot2)
##use ggpattern to get rid of black color in graph
library(ggpattern)

##transform data
transformed_data <- transform(
  top10_melt,
  X2 = factor(X2, levels = c("HepaRG", "RPTEC", "combined")),
  Categories = factor(Categories, levels = row.names(top10)),
  pattern = ifelse(
    factor(top10_melt$X3, levels = rev(c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "2Phenylphenol"))) %in% 
      rev(c("Chlorotoluron", "Thiabendazole", "2Phenylphenol")), "stripe", "none"
  ))

## Create the plot
gg <- ggplot(
  transformed_data,
  aes(
    y = -log10(value),
    x = Categories,
    fill = factor(X3, levels = rev(c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "2Phenylphenol"))),
    pattern = pattern
  )
)

## Add geom_col_pattern
gg <- gg + 
  geom_col_pattern(
    width = 0.6,
    position = position_dodge(width = 0.7),
    color = "black",
    lwd = 0.8,
    pattern_colour = "grey34",
    pattern_fill = "grey34",
    pattern_density = 0.4,  # Adjust stripe density
    pattern_spacing = 0.02  # Adjust stripe spacing
  )+
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"))

## Add coordinates and facets
gg <- gg +
  coord_flip() +
  scale_x_discrete(labels = gsub(",", "\n", row.names(top10))) +
  facet_wrap(~X2, ncol = 3)

## Add guides and scales
gg <- gg +
  scale_fill_manual(
    values = c("white", "grey74", "grey24", "white", "grey74", "black"),
    breaks = c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "2Phenylphenol"),
    labels = c("Cyproconazole", "Fluxapyroxad", "Azoxystrobin", "Chlorotoluron", "Thiabendazole", "2-Phenylphenol")
  ) 

## Add labels and theme
gg <- gg +
  labs(x = "", y = expression("-"*log[10]*"(p-value)")) +
  guides(
    pattern = "none",
    fill = guide_legend(
      title = "",
      override.aes = list(
        pattern = c("none", "none", "none", "stripe", "stripe", "stripe"),
        fill = c("white", "grey74", "grey24", "white", "grey74", "black")
      )
    )
  )+
  theme_bw()

ggsave(
  filename = "IPA.png", device = "png",
  width = 297,
  height = 210,
  units = "mm",
  dpi = 900)
