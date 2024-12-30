library(dplyr)
library(tidyr)
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
setwd("C:/Users/ASUS/OneDrive/Desktop/Project/Tomato_r")
count_read<- read.csv("GSE220689_Transcript_read_count_all_conditions.csv")

count_tab<- count_read%>%
  subset(., select= c(2,4,5,6,7,8,9,10,11))


# get metadata
gse<- getGEO(GEO = "GSE220689", GSEMatrix = TRUE)
metadata<- pData(phenoData(gse[[1]]))
head(metadata)

metadata.mod<- metadata%>%
  subset(., select=c(1,37,39))


metadata.mod<- metadata.mod%>%
  mutate(title= gsub("2 hours ", "T2_", title))%>%
  mutate(title= gsub("6 hours ", "T6_", title))%>%
  mutate(title=gsub("control ", "control", title))%>%
  mutate(title= gsub("tunicamycin ", "tunik", title))%>%
  mutate(title= gsub("biol ", "", title))%>%
  mutate(title= gsub("rep 1", "_rep1", title))%>%
  mutate(title= gsub("rep 2", "_rep2", title))
colnames(metadata.mod)

#save the count table
write.csv(count_tab, "raw_counts.csv")

typeof(metadata.mod$title)
 design<- read.csv("design.csv", row.names = 1)
 unique(design$Group)

# now we do according  to Bioinformatics coach
 





