library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(glmGamPoi)

mydata = SCTransform(mydata, method="glmGamPoi", vars.to.regress="percent.mito", verbose=FALSE)
mydata = RunPCA(mydata, verbose=FALSE)
mydata = RunHarmony(mydata, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(mydata)
mydata <- FindNeighbors(mydata, dims=1:20, reduction="harmony")
mydata <- RunUMAP(mydata, dims=1:20, reduction="harmony")

colors = c("#66CCCC", "#FF99CC", "#CCFF66")

######################################
DimPlot(mydata, reduction="umap", group.by="patient_ID", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+
scale_color_manual(values=colors)
#################################################################################
library(clustree)
obj = FindClusters(mydata, resolution = seq(0.1, 0.3,by=0.05))
clustree(obj)
#############################
mydata <- FindClusters(mydata, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
#################################################################################
FeaturePlot(mydata, features=c("Ly6d"), cols=c("lightgray", "red"))+NoLegend()
VlnPlot(mydata, features=c("AGER"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0: Monocytes: OLR1, AQP9
# 1: M2 Macrophages: CCL18, FOLR2, CD163
# 2: cDCs 1: CLEC9A, WDFY4
# 3: Basal cells: DSC3, KRT6B, LGALS7, DSG1
# 4: cDCs 2: CCL22, CCL17, CD1B
# 5: Fibroblasts: DCN, COL1A1, LUM
# 6: Langerhans cells: CD207, CD1A
cell_label = c("Monocytes", "M2 Macrophages", "cDCs 1", "Basal cells", "cDCs 2", "Fibroblasts", "Langerhans cells")
#################################################################################
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)

set = c("OLR1", "AQP9", "CCL18", "FOLR2", "CD163", "CLEC9A", "WDFY4", "DSC3", "KRT6B", "LGALS7", "CCL22", "CCL17", "CD1B", "DCN", "LUM", "CD207", "CD1A")
DotPlot(mydata, features=set)+coord_flip()+scale_color_distiller(palette="RdYlBu")+theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())

###########################################################################
bar = mydata@meta.data %>% group_by(Type, cell_type) %>% count()
bar$cell_type = factor(bar$cell_type, levels=cell_label)
bar$Type = factor(bar$Type, levels=c("Normal", "Tumor", "Liver_metastasis"))
ggplot(data=bar, aes(x=Type, y=n, fill=cell_type))+ 
geom_bar(stat="identity", position=position_fill())+
scale_fill_manual(values=colors)+theme_classic()+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=30, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")

#####################################################################
Tissue_label = c("primary tumor", "metastasis")
bar2$Tissue = factor(bar2$Tissue, levels=Tissue_label)
ggplot(data=bar2, aes(x=cell_type, y=percent, fill=Tissue))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c("dodgerblue1", "darkorange1"))+theme_classic()+geom_text(aes(label=percent), vjust=-0.2)+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())

############################################################################
MMP12 = subset(mydata, cell_type=="MMP12+ mydata")
Idents(MMP12) = MMP12@meta.data$Type
DEG <- FindAllMarkers(MMP12, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(DEG, "MMP12_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
saveRDS(mydata, "./mydata_sub.rds")






