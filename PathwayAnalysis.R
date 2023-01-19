BiocManager::install("clusterProfiler")
BiocManager::install("Affyhgu133A2Expr")
library(clusterProfiler)
library(tidyverse)
data("gcSample")

setwd("/Users/koji11235/Project/RNAseqAnalysis_airway/")
read_csv("significant_genes_airway.csv") %>% 
  filter(padj<1e-2 & abs(log2FoldChange)>2) %>% 
  write_csv("significant_genes.csv")

sample_data<-read.csv("significant_genes.csv")
differential_expressed_genes<-sample_data$Entrez_Gene_ID
kegg_enrichment_result <- enrichKEGG(gene=differential_expressed_genes,pvalueCutoff=0.05)
as.data.frame(kegg_enrichment_result)->a


#ライブラリの読み込み
library(clusterProfiler)

#サンプルデータの読み込み
sample_data<-read.csv("sample_gene_data.csv")
#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID
#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)
#結果の表示
head(as.data.frame(kegg_enrichment_result))
head(kegg_enrichment_result, n=10)

browseKEGG(kegg_enrichment_result,"hsa04510")

#BiocManager::install("pathview")
#ライブラリの読み込み
library(pathview)

#サンプルデータからlog2FoldChangeを抽出
logFC <- sample_data$log2FoldChange

#logFCとEntrez_Gene_IDを対応させる
names(logFC) <- sample_data$Entrez_Gene_ID

#パスウェイの画像のダウンロード
pathview(gene.data = logFC, pathway.id = "hsa04510", limit = list(gene=2, cpd=1))



BiocManager::install("ReactomePA")
#ライブラリの読み込み
library(ReactomePA)

#サンプルデータの読み込み
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

#サンプルデータからlog2FoldChangeを抽出
logFC <- sample_data$log2FoldChange

#logFCとEntrez_Gene_IDを対応させる
names(logFC) <- sample_data$Entrez_Gene_ID

#Reactomeパスウェイ解析を実行
Reactome_enrichment_result <- enrichPathway(gene=gene_IDs,pvalueCutoff=0.05, readable=T)

#結果の表示
head(as.data.frame(Reactome_enrichment_result))
#バープロット
barplot(Reactome_enrichment_result, showCategory=8,x = "Count")
#ドットマッププロット
dotplot(Reactome_enrichment_result, showCategory=15)
#エンリッチメントマッププロット
emapplot(Reactome_enrichment_result)
#cnetプロット
cnetplot(Reactome_enrichment_result, categorySize="pvalue", foldChange=logFC,showCategory = 8)



emapplot(kegg_enrichment_result)


#goplot(Reactome_enrichment_result)
cnet<-cnetplot(Reactome_enrichment_result, categorySize="pvalue", foldChange=logFC,showCategory = 8)
class(cnet)
bar<-barplot(Reactome_enrichment_result, showCategory=8,x = "Count")
class(bar)
ggsave("bartest.png")

#パスウェイIDを取得（R-HSA-9006934）
pathwayID<-Reactome_enrichment_result[2]$ID

#URLを設定
URL<-paste0("https://reactome.org/ContentService/exporter/diagram/",pathwayID,".png")

#保存するファイル名を設定
output_filename<-paste0(pathwayID,".png")

#ダウンロードを実行
download.file(URL,output_filename)

#文字列加工用のライブラリー
library(stringr)

#パスウェイIDを取得（R-HSA-9006934）
pathwayID<-Reactome_enrichment_result[2]$ID

#含まれている遺伝子を抽出，加工
genelist<-str_split(Reactome_enrichment_result[2]$geneID,"/")[[1]]
genelist.str<-str_c(genelist,collapse = ",")

#URLを設定
URL<-paste0("https://reactome.org/ContentService/exporter/diagram/",pathwayID,".png","?flg=",genelist.str)

#保存するファイル名を設定
output_filename<-paste0(pathwayID,"_decorated",".png")

#ダウンロードを実行
download.file(URL,output_filename)

#############################################################
#パスウェイIDを取得（R-HSA-9006934）
pathwayID<-Reactome_enrichment_result[2]$ID

genelist<-str_split(Reactome_enrichment_result[2]$geneID,"/")[[1]]
Entrez_Gene_ID<-mapIds(org.Hs.eg.db, keys=genelist, column='UNIPROT', keytype='SYMBOL')
Entrez_Gene_ID_list<-str_c(Entrez_Gene_ID,collapse = ",")

#URLを設定
URL<-paste0("https://reactome.org/ContentService/exporter/diagram/",pathwayID,".png","?flg=",Entrez_Gene_ID_list)
URL<-paste0("https://reactome.org/ContentService/exporter/diagram/",pathwayID,".png","?flg=",genelist %>% str_c(collapse = ","))

#保存するファイル名を設定
output_filename<-paste0(pathwayID,"_decorated_2",".png")

#ダウンロードを実行
download.file(URL,output_filename)


i=2

#パスウェイIDを取得（R-HSA-9006934）
pathwayID<-Reactome_enrichment_result[i]$ID

genelist<-str_split(Reactome_enrichment_result[i]$geneID,"/")[[1]]
Entrez_Gene_ID<-mapIds(org.Hs.eg.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
Entrez_Gene_ID_list<-str_c(Entrez_Gene_ID,collapse = ",")

#URLを設定
URL<-paste0("https://reactome.org/ContentService/exporter/diagram/",pathwayID,".png","?flg=",Entrez_Gene_ID_list)

#保存するファイル名を設定
output_filename<-paste0(pathwayID,"_decorated_",".png")

#ダウンロードを実行
download.file(URL,output_filename)

