library(ggplot2)
qplot(
  unlist(x[x$gene_name_column == "A1BG", -1]), main = "Expression of A1BG", xlab = "Expression", ylab = "Count")

#
a <- as.numeric(unlist(aggregated_expr[aggregated_expr$Gene == "AAAS", -1]))
b <- as.numeric(unlist(aggregated_expr[aggregated_expr$Gene == "7F5", -1]))


#qplot
qplot(unlist(aggregated_expr["AAAS", -1]))

# qplot(a)
qplot(a,
      b,
      colour = ("red"),
      xlab = "AAAS",
      main = "Expression of AAAS")

aggregated_expr
data_cl <- matrix( sample(seq(1,2000),200), ncol = 10 )
rownames(data) <- paste0("sample_" , seq(1,20))
colnames(data) <- paste0("variable",seq(1,10))

clust1 <- na.omit(aggregated_expr[1:20,])
clust1 <- scale(clust1)
dist(clust1 , method ="euclidean")
dist_1 <- dist(t(df_for_heatmap), c(4:8) , diag=TRUE)
hc <- hclust(dist_1)
plot (hc)

hc1 <- hclust(df_for_heatmap, method = "complete" )



