#https://programming-workshops.readthedocs.io/en/stable/workshops/04_R/Workshop_R_Solution.html
setwd("/Users/koji11235/Project/RNAseqAnalysis_airway/")
library(edgeR)

BiocManager::install(c("airway","DESeq2"))
library(tidyverse)
library(magrittr)
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(airway)

data("airway")

head(assay(airway))
DE_airway <- DESeqDataSet(airway, design = ~ cell + dex)
DE_airway
DE_airway@colData$dex<-relevel(DE_airway@colData$dex, ref = "untrt")

help("DESeq")

DE_airway <- DESeq(DE_airway)
res <- results(DE_airway)
res
res[order(res$pvalue),]
plotMA(res, ylim=c(-2,2))

res@listData
res.filtered<-as.data.frame(res) %>% 
  rownames_to_column("Ensembl") %>% 
  as_tibble() %>% 
  mutate(sig=padj<0.01) %>% 
  filter(padj<0.01) 
res.filtered%>%
  #mutate(ENTREZID=clusterProfiler::bitr(res.filtered$Ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = FALSE) %>% .$ENTREZID )
  mutate(Entrez_Gene_ID=mapIds(org.Hs.eg.db, keys=Ensembl, column='ENTREZID', keytype='ENSEMBL')) %>% 
  dplyr::select(Ensembl,Entrez_Gene_ID,everything(),-sig) %>% 
  write_csv("significant_genes_airway.csv")

as.data.frame(res) %>% 
  rownames_to_column("Ensembl") %>% 
  as_tibble() %>% 
  mutate(sig=padj<0.01) %>% 
  ggplot(aes(log2FoldChange, -1*log10(padj), col=sig)) + geom_point() + ggtitle("Volcano plot")

as.data.frame(res) %>% 
  rownames_to_column("Ensembl") %>% 
  as_tibble() %>% 
  mutate(sig=padj<0.01) %>% 
  ggplot( aes(baseMean, log2FoldChange, col=sig)) + geom_point() + scale_x_log10() + ggtitle("MA plot")






