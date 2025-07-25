devtools::install_github("sqjin/CellChat")
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(svglite)
options(stringsAsFactors=FALSE)

mydata = readRDS("./Neuron3.rds")
cellchat = createCellChat(object=mydata, group.by="cell_type", meta=mydata@meta.data)
groupSize = as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.human

# > unique(CellChatDB$interaction$annotation)
# [1] "Secreted Signaling"	"ECM-Receptor"	"Cell-Cell Contact"

CellChatDB.use <- subsetDB(CellChatDB, search="Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = projectData(cellchat, PPI.human)

#####################################################################################
cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)
cellchat = filterCommunication(cellchat, min.cells=10)
df.net = subsetCommunication(cellchat)
cellchat = computeCommunProbPathway(cellchat)
df.netp = subsetCommunication(cellchat, slot.name="netP")

#####################################################################################
cellchat = aggregateNet(cellchat)
groupSize = as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight=groupSize, weight.scale=T, label.edge=T, title.name="Number of interactions", arrow.size=1, color.use=c("hotpink", "deepskyblue", "forestgreen"), edge.label.cex=1.2)
netVisual_circle(cellchat@net$weight, vertex.weight=groupSize, weight.scale=T, label.edge=T, title.name="Number of strength", arrow.size=1, color.use=c("hotpink", "deepskyblue", "forestgreen"))

mat = cellchat@net$count
par(mfrow=c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)){
	mat2 = matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
	mat2[i, ] = mat[i, ]
	netVisual_circle(mat2, vertex.weight=groupSize, weight.scale=T, arrow.width=0.2, arrow.size=0.1, edge.weight.max=max(mat), title.name=rownames(mat)[i])
}

#########################################################################################
netAnalysis_contribution(cellchat, signaling=pathways.show)
pairLR.TGFb = extractEnrichedLR(cellchat, signaling=pathways.show, geneLR.return=FALSE)
#########################################################################################
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use=c(3,5,7,8,9), targets.use=c(1,2,4,6), remove.isolate=FALSE)
netVisual_bubble(cellchat, sources.use=c(3,5,7,8,9), targets.use=c(1,2,4,6), signaling=c("CCL", "CXCL"), remove.isolate=FALSE)

netVisual_bubble(cellchat, signaling=c("IL2"), remove.isolate=FALSE)+theme(axis.text.x=element_text(angle=15, hjust=1))

#################################################################################################
#################################################################################################
#################################################################################################
cellchat_Baseline = createCellChat(object=mydata_Baseline, group.by="cell_type", meta=mydata_Baseline@meta.data)
cellchat_Day14 = createCellChat(object=mydata_Day14, group.by="cell_type", meta=mydata_Day14@meta.data)

cellchat_Baseline@DB = CellChatDB.human
cellchat_Baseline = subsetData(cellchat_Baseline)
cellchat_Baseline = identifyOverExpressedGenes(cellchat_Baseline)
cellchat_Baseline = identifyOverExpressedInteractions(cellchat_Baseline)
cellchat_Baseline = projectData(cellchat_Baseline, PPI.human)
cellchat_Baseline = computeCommunProb(cellchat_Baseline, raw.use=FALSE, population.size=TRUE)
cellchat_Baseline = filterCommunication(cellchat_Baseline, min.cells=10)
cellchat_Baseline = computeCommunProbPathway(cellchat_Baseline)
cellchat_Baseline = aggregateNet(cellchat_Baseline)

cellchat_Day14@DB = CellChatDB.human
cellchat_Day14 = subsetData(cellchat_Day14)
cellchat_Day14 = identifyOverExpressedGenes(cellchat_Day14)
cellchat_Day14 = identifyOverExpressedInteractions(cellchat_Day14)
cellchat_Day14 = projectData(cellchat_Day14, PPI.human)
cellchat_Day14 = computeCommunProb(cellchat_Day14, raw.use=FALSE, population.size=TRUE)
cellchat_Day14 = filterCommunication(cellchat_Day14, min.cells=10)
cellchat_Day14 = computeCommunProbPathway(cellchat_Day14)
cellchat_Day14 = aggregateNet(cellchat_Day14)
cellchat_list = list(Baseline=cellchat_Baseline, Day14=cellchat_Day14)
cellchat_combine = mergeCellChat(cellchat_list, add.names=names(cellchat_list), cell.prefix=T)
compareInteractions(cellchat_combine, show.legend=T, group=c(1, 2), measure="count", color.use=c("#CCFF66", "#FF99CC"))
rankNet(cellchat_combine, mode="comparison", stacked=T, do.stat=T, color.use=c("#CCFF66", "#FF99CC"))



