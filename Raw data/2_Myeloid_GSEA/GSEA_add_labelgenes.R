devtools::install_github("junjunlab/GseaVis")
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GseaVis)

mydata = readRDS("../1_preprocessing/mydata_cluster.rds")
mydata = subset(mydata, cell_type=="Myeloid cells")
data = FindMarkers(mydata, group.by="Type", ident.1="metastasis", ident.2="primary tumor", logfc.threshold=0, min.pct=0)
write.table(data, "./Myeloid_treatment_control_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")

symbol <- rownames(data)
entrez <- bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(entrez)
#    SYMBOL ENTREZID
# 11   A1BG        1
# 12   A1CF    29974
# 13    A2M        2
# 14 A4GALT    53947
# 15  A4GNT    51146
# 16   AAAS     8086

genelist <- data$log2FoldChange
names(genelist) <- rownames(data)
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist), entrez[,1]), 2]
genelist = sort(genelist, decreasing=T)

##############################################################################################
R.utils::setOption("clusterProfiler.download.method", "auto")
KEGG_ges <- gseGO(
geneList=genelist,
ont="BP",
OrgDb=org.Hs.eg.db,
keyType="SYMBOL",
minGSSize=10,
maxGSSize=500,
pvalueCutoff=0.05,
pAdjustMethod="BH",
verbose=FALSE,
eps=0
)
write.table(KEGG_ges@result, "./GOBP_result.txt", col.names=T, row.names=T, quote=F, sep="\t")
#################################################################################
KEGG_ges2 <- setReadable(KEGG_ges, OrgDb=org.Hs.eg.db, keyType="SYMBOL")
core <- result$core_enrichment[58]
core_genes <- str_split(core ,'/')[[1]]
p = gseaNb(
    object = KEGG_ges,
    geneSetID = result$ID[58],
    addPval = T,
    pvalX = 0.95,
    pvalY = 0.8,
    newGsea = T,
    addPoint = F,
    newCurveCol = c("steelblue1","lightgray", "hotpink1"),
    newHtCol = c("hotpink1", "lightgray", "steelblue1"),
    addGene = core_genes, 
    kegg = F,
    geneCol = 'black'
)
ggsave("PGSEA_Keratinization.pdf", p, width=6, height=6)
##########################################################################
library(ggplot2)
df = read.table("./GSEA_GOBP_result_barplot.txt", header=T, sep="\t")
ggplot(df, aes(x=NES, y=reorder(Description, NES), color=pvalue))+
geom_point(aes(size=abs(NES)))+
geom_segment(aes(x=0, xend=NES, y=Description, yend=Description), color="black", linewidth=1)+
scale_color_gradient(low="hotpink1", high="steelblue1")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=15))




