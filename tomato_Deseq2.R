rm(list = ls())

library(DESeq2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
setwd("C:/Users/ASUS/OneDrive/Desktop/Project/Tomato_r")

count_table<- read.csv("raw_counts.csv")
dim(count_table)

#having problem with duplicates
# Check for duplicate row names 
duplicated_rows <- count_table[duplicated(count_table[,1]), ] 
print(duplicated_rows) 
# Remove duplicates
count_table <- count_table[!duplicated(count_table[,1]), ]
# Set row names
rownames(count_table) <- count_table[,1] 
count_table <- count_table[,-1] 
# View the updated data frame 
print(head(count_table))

sample_info<- read.csv("design.csv", row.names = 1)
dim(sample_info)

#making sure that row names in sample_info matches colnames in count_table
all(colnames(count_table)%in% rownames(colData))

colnames(count_table)
rownames(sample_info)

# Check for exact matches
all(trimws(colnames(count_table)) %in% trimws(rownames(sample_info)))

# remove extra space
rownames(sample_info)<- trimws(rownames(sample_info))
colnames(count_table)<- trimws(colnames(count_table))
colnames(sample_info)<- trimws(colnames(sample_info))
all(colnames(count_table)%in% rownames(sample_info))


# setting factor level
sample_info$Group<- as.factor(sample_info$Group)
sample_info$Time<- as.factor(sample_info$Time)
groups<- unique(sample_info$Group)

print(colnames(sample_info))
colnames(sample_info) <- trimws(colnames(sample_info))
print(unique(sample_info$Group))

# create DESeq object
dds <- DESeqDataSetFromMatrix( countData = count_table, colData = sample_info, design = ~ Group )
print(dds)

#pre-filtering
keep<- rowSums(counts(dds))>=10
dds<- dds[keep,]

print(dds)

# setting factor level

dds$Group<- relevel(dds$Group, ref = "control")
# note : collapse technical replicates


# running DESeq2
dds<- DESeq(dds)
res<- results(dds)
res

# changing DESeq object to R object (data frame)

res_data<- as.data.frame(res)

#explore data frame
summary(res_data)

#order the data frame to p-adj
# NA value represents gene expression could not determined by "DESeq2"
res_data_order<- res_data[order(res_data$pvalue),]
head(res_data_order)

# explore results
summary(res)

#adjusted p value threshold of 0.01
res0.01<- results(dds, alpha = 0.01)
summary(res0.01)

#contrasts
resultsNames(dds)

#Visualize the plot
#MA plot
plotMA(res)


# is ZIP differentially expressed?
res_data["ZIP",]

# selecting genes with a significant change in gene expression (adjusted p-value below 0.05)
# and log2FoldChange <1 and >1



#step 1
filtered<- res_data%>%
  filter(res_data$padj< 0.05)
# step 2
filtered<- filtered%>%
  filter(abs(filtered$log2FoldChange)>1)
# checking how much we filtered
dim(res_data)
dim(filtered)


# saving the results
write.csv(res_data, "result_data_table.csv")
write.csv(filtered, "filter_result_data.csv")

# normalize the read counts
normalized_counts<- counts(dds, normalized= TRUE)
write.csv(normalized_counts, "normalized_counts.csv")

# visualization
plotDispEsts(dds)

# Estimate dispersions 
dis <- estimateDispersions(dds)
# Plot dispersions 
plotDispEsts(dis)


#PCA Plots

# variance stabilizing transformation
vsd<- vst(dds, blind = FALSE)

# use this transformed valu for PCA
plotPCA(vsd, intgroup= "Group")


#pheat map library
library(pheatmap)
# heatmap  of sample-to-sample distance matrix( with clustering) based on the normalized counts

#generate distance matrix
sampleDists<- dist(t(assay(vsd)))
sampleDistsMatrix<-  as.matrix(sampleDists)
colnames(sampleDistsMatrix)


# heatmap
pheatmap(sampleDistsMatrix, main = "Heat Map of Sample Distance Matrix Clustering", 
         clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)



# heat map of log transformed normalized counts, using top 10 genes

#top 10 genes
top_hits<- res_data[order(res_data$padj),][1:10,]
head(top_hits)
top_hits<- row.names(top_hits)
top_hits

rld<-rlog(dds, blind = FALSE)

# heat maps with out clustering
pheatmap(assay(rld)[top_hits,], cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE)

# add details to the graph
# annot_info<- as.data.frame(colData(dds)[,c("Group")])

# heat map with clustering
pheatmap(assay(rld)[top_hits,], main = "Heat Map of Top 10 log Transformed Normalized Counts of Genes in Different Samples")
  


# heat map of Z scores

# claculate Z scors
cal_Zscore<- function(x){(x-mean(x)/ sd(x))}

zscore_all<- t(apply(normalized_counts, 1,cal_Zscore))
zscore_subset<- zscore_all[top_hits,]
pheatmap(zscore_subset, main="Heat Map of Top 10 Z-Score Normalized Counts of Genes in Different Samples")
  


# MA Plot
plotMA(dds, ylim= c(-3,3))

# remove noise
resLFC<- lfcShrink(dds, coef = "Group_tunik_vs_control", type = "apeglm")
plotMA(resLFC, title="MA Plot with Reduced Noise Using Apeglm")



# volcano plot

#change resLFC to data frame

resLFC<- as.data.frame(resLFC)
install.packages("ggrepel")
library(ggrepel)

# label the genes

resLFC$diffexpressed<-"NO"
resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05]<- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05]<- "Down"
resLFC$delabel<- NA

ggplot(dat= resLFC, aes(x= log2FoldChange,y= -log10(pvalue), col= diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c("blue","black", "red"))+
  theme(text = element_text(size = 20))
























































