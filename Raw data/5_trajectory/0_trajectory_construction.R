library(monocle)
mydata = readRDS("./Epithelial_subcluster.rds") 
expr_matrix = as(as.matrix(mydata$RNA@counts), "sparseMatrix")
p_data = mydata@meta.data
f_data = data.frame(gene_short_name=row.names(expr_matrix), row.names=row.names(expr_matrix))
pd = new("AnnotatedDataFrame", data=p_data)
fd = new("AnnotatedDataFrame", data=f_data)

cds = newCellDataSet(expr_matrix, phenoData=pd, featureData=fd, lowerDetectionLimit=0.5, expressionFamily=negbinomial.size())

######################################################################################
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

cds = detectGenes(cds, min_expr=0.1)
expressed_genes = row.names(subset(fData(cds), num_cells_expressed>=10)) 

# step 1: choosing genes that define progress
# step 2: reducing the dimensionality of the data
# step 3: ordering the cells in pseudotime

# disp_table = dispersionTable(cds)
# disp.genes = subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# cds = setOrderingFilter(cds, disp.genes)
########################################################################

expressed_genes = VariableFeatures(mydata)
# diff = differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr="~Type", cores=1)
qval<0.01,decreasing=F

Idents(mydata) = mydata@meta.data$Type
DEG = FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
DEG = subset(DEG, p_val_adj<0.01)
DEG = subset(DEG, avg_log2FC>0.25)
ordergene = unique(DEG$gene)
cds = setOrderingFilter(cds, ordergene)
# plot_ordering_genes(cds)
cds = reduceDimension(cds, max_components=2, method="DDRTree")
cds = orderCells(cds, root_state=5)

########################################################################
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=T)

colors = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
p = plot_cell_trajectory(cds, color_by="Sample", size=1, show_backbone=T)+scale_color_manual(values=colors)
ggsave("./trajectory_Type.pdf", p, width=6, height=6)

p = plot_cell_trajectory(cds, color_by="Pseudotime")+scale_color_gradient(low="lightgrey", high="red")
ggsave("./trajectory_pseudotime.pdf", p, width=6, height=6)

#######################################################################################
library(ggpubr)
df <- pData(cds) 

p = ggplot(df, aes(x=Pseudotime, color=Type, fill=Type))+
scale_color_manual(values=c("dodgerblue", "hotpink"))+
scale_fill_manual(values=c("dodgerblue", "hotpink"))+
geom_density(alpha=0.8)+
theme_classic()+
theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.x=element_text(size=14, face="bold"), axis.title.y=element_text(size=14, face="bold"), legend.text=element_text(size=14, face="bold"), legend.title=element_blank(), legend.position="top")
ggsave("./density_Type.pdf", p, height=3, width=6)

set = c("RELN", "MSRA", "SYN2", "PSEN1")
plot_genes_in_pseudotime(cds[set, ], color_by="Type", cell_size=1)+scale_color_manual(values=c("dodgerblue", "hotpink"))

genes = c("PSEN1", "DAAM2", "DCC", "CLU")
plot_genes_branched_pseudotime(cds[genes, ], branch_point=1, color_by="State", cell_size=1, ncol=1, branch_labels=c("branch 1", "branch 2"))+scale_color_manual(values=c("orange", "lightgray", "forestgreen"))

##########################################################
ggplot(df, aes(x=Pseudotime, y=GOBP_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY1))+
geom_point(alpha=1, size=2, color="orange")+
geom_smooth(method="lm", color="black")+
theme_bw()+
theme(axis.text.x=element_text(face="bold", size=15), axis.text.y=element_text(face="bold", size=15), axis.title=element_text(face="bold", size=15), legend.position="none")
cor.test(df$Pseudotime, df$GOBP_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY1, alternative="two.side", method="pearson", conf.level=0.95)
# apoptotic:      R=0.3311, P<2.2e-16
# neuron death:   R=-0.0513, P=0.02961
# WNT pathway:    R=0.5369, P<2.2e-16
############################################################################################################
genes = c("HLA-A", "CD151", "PAK1", "EPCAM")
colors = c("chocolate1", "peachpuff1", "olivedrab1", "olivedrab1", "olivedrab1")
plot_genes_branched_pseudotime(cds[genes, ], branch_point=2, color_by="State", cell_size=1, ncol=1, branch_labels=c("Cell fate 1", "Cell fate 2"))+scale_color_manual(values=colors)

set = c("TCF7L2", "WNT2B", "CCNY", "CSNK1A1", "CTNNBIP1", "MCC", "CDK14")
plot_genes_in_pseudotime(cds[set, ], color_by="Diagnosis", cell_size=1)+scale_color_manual(values=c("#9EB9F3", "#FC8D62"))



