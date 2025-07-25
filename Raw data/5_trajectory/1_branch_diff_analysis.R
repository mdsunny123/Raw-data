# (1)Basic Differential Analysis
# (2)Finding Genes that Distinguish Cell Type or State
# (3)Finding Genes that Changes as a Function of Pseudotime
ordergenes = as.vector(disp.genes)
Time_diff = differentialGeneTest(cds[ordergenes, ], cores=1, fullModelFormulaStr="~sm.ns(Pseudotime)")
Time_diff = subset(Time_diff, qval<0.01)
Time_diff = arrange(Time_diff, qval)
Time_diff = Time_diff[, c(5,2,3,4,1,6,7)]
Time_genes = Time_diff %>% pull(gene_short_name) %>% as.character()

plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters=3, show_rownames=T, return_heatmap=T)

p = plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters=3, show_rownames=F, return_heatmap=T)
ggsave("./heatmap_pseudotime.pdf", p, width=4.5, height=6)
clusters = cutree(p$tree_row, k=2)
clustering = data.frame(clusters)
clustering[, 1] = as.character(clustering[, 1])
colnames(clustering) = "Gene_Clusters" 
write.table(clustering, "./pseudotime_clusters.txt", col.names=T, row.names=T, sep="\t", quote=F)

######################################################################################################
BEAM_res = BEAM(cds[ordergene, ], branch_point=2, cores=1)
BEAM_res = BEAM_res[order(BEAM_res$qval), ]
BEAM_res = BEAM_res[, c("gene_short_name", "pval", "qval")]
BEAM_res = subset(BEAM_res, qval<0.01)
p = plot_genes_branched_heatmap(cds[row.names(BEAM_res),], branch_point=2, num_clusters=3, show_rownames=F, branch_labels=c("Cell fate 1", "Cell fate 2"), return_heatmap=T, branch_colors=c("peachpuff1", "chocolate1", "olivedrab1"))
###################################################
branch_clusters = cutree(p$ph_res$tree_row, k=3)
branch_clustering = data.frame(branch_clusters)
branch_clustering[, 1] = as.character(branch_clustering[, 1])
colnames(branch_clustering) = "Branch_Clusters"
write.table(branch_clustering, "./BEAM_clusters.txt", col.names=T, row.names=T, sep="\t", quote=F)


set = c("PGK1", "DDX3Y", "JUND", "HSPD1", "FTH1")
plot_genes_in_pseudotime(cds[set, ], color_by="Type", cell_size=1)+scale_color_manual(values=c("dodgerblue", "hotpink"))


genes = c("MFAP4", "VCAN", "CXCL12", "ITGA10", "PCDH7", "FN1", "THBS2", "MYH10", "ENG")
genes = c("CDKN1A", "BTG2", "BCL6", "ADAMTS1", "HSPA1B", "HSPA1A", "FHL1", "JADE1")

