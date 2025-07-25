library(ggplot2)
df = read.table("./GSEA_GOBP_result_barplot.txt", header=T, sep="\t")
p = ggplot(df, aes(x=NES, y=reorder(Description, NES), fill=pvalue))+
geom_bar(stat='identity', width=0.6)+
scale_fill_gradient(low="#FF8C00", high="#00BFFF")+
theme_bw()+
theme(axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=15))
ggsave("./GSEA_GOBP_result_barplot.pdf", p, width=6, height=6)



