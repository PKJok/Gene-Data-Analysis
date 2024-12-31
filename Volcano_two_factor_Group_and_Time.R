# volcano plot of two factor "group" and "time

rm(list = ls())

library(DESeq2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
setwd("C:/Users/ASUS/OneDrive/Desktop/Project/Tomato_r")
reslfc2<- read.table("resLFC_two_factor.csv", sep=",", header= TRUE)
head(reslfc2)

colnames(reslfc2)[1]<- "Gene_Name"
head(reslfc)

#set Theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

reslfc2$diffexpressed<-"NO"
reslfc2$diffexpressed[reslfc2$log2FoldChange < -0.6 & reslfc2$pvalue < 0.05]<-"Down"
reslfc2$diffexpressed[reslfc2$log2FoldChange > 0.6 & reslfc2$pvalue < 0.05]<-"UP"

reslfc2$delabel <- ifelse(reslfc2$Gene_Name %in% head(reslfc2[order(reslfc2$padj), "Gene_Name"], 20),
                         reslfc2$Gene_Name, NA)

# making volcano plot
ggplot(reslfc2, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label= delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#4682B4", "#778899", "maroon"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 20), xlim = c(-3.5, 3.5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Treatment (control vs tunicamycin) Volcano Plot')+ # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 
