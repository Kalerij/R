#to Russian language
Sys.setlocale('LC_ALL', 'russian')
####
# setwd("C:/Users/rodic/OneDrive/Документы/ ")
gene_set_hgnc1<- read.delim("C:/Users/rodic/OneDrive/Документы/GPL14951-11332.txt",fill = T,skip=28)

platform<- read.table("C:/Users/rodic/OneDrive/Документы/GPL14951-11332.txt",fill = T,quote = "\"",header=T,skip = 28,sep="\t")

expression_probes<- read.csv("C:/Users/rodic/OneDrive/Документы/GSE52219_series_matrix.txt",fill = T, skip=62,sep="\t")

expression_probes$ILMN_GENE <-
  platform$ILMN_Gene[match(expression_probes$ID_REF, platform$ID)]
length(expression_probes$ILMN_GENE)
for(i in 2:(ncol(expression_probes)-1)){
  expression_probes[,i]<-as.numeric(expression_probes[,i])
}
expression_probes[is.na(expression_probes)]<-0
sum(is.na(expression_probes))

library(psych)

aggregated_expr<-aggregate(expression_probes[,-c(1,ncol(expression_probes))], by=list (expression_probes$ILMN_GENE),FUN = "geometric.mean")

head(aggregated_expr)

colnames(aggregated_expr)[1]<-"Gene"
head(aggregated_expr)
rownames(aggregated_expr)<-aggregated_expr$Gene

library(ggplot2)

number_of_samples<-ncol(aggregated_expr)

df_long_for_2_genes<- data.frame("gene_name_column"=c(rep("A2M",number_of_samples),rep("A2LD1",number_of_samples)),"expression"=c(unlist(aggregated_expr["A2M",]),unlist(aggregated_expr["A2LD1",])))
qplot(unlist(x[x$gene_name_column == "A2M",-1]), main = "Expression of A2M", xlab = "Expression", ylab = "Count")

plot(
  x = unlist(aggregated_expr[aggregated_expr$Gene == "AAAS", -1]),
  y = unlist(aggregated_expr[aggregated_expr$Gene == "A2M", -1]),
  main = "expression for AAA1 vs A2M",
  ylab = "AAAS",
  xlab = "A2M",
  cex = 1,
  col = "red",
  pch = 19
)
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

# 31.10
aggregated_expr<-aggregated_expr[,-1]
head(aggregated_expr)
ncol(aggregated_expr)
group1 <-1:12
group2 <-13:23
  

#   
t.test(unlist(aggregated_expr[2,group1]),unlist(aggregated_expr[2,group2]))

wilcox.test(unlist(aggregated_expr[2,group1]),unlist(aggregated_expr[2,group2]))

# можно
pv <- NA [1:10]
pv_w <- NA[1:10]
library(forecast)
for (i in 1:100){
  gene_expr_for_resp<-unlist(aggregated_expr[i,group1])
  gene_expr_for_nonresp<-unlist(aggregated_expr[i,group2])

  # if both vectors is not constants and not na
  
  if(!is.constant(gene_expr_for_resp) & !is.constant(gene_expr_for_nonresp)& sum(is.na(gene_expr_for_nonresp))==0 & sum(is.na(gene_expr_for_resp))==0 ){
    print(i)
  t1<-
  t.test(gene_expr_for_resp,gene_expr_for_nonresp)
pv[i] <-t1$p.value
t2<-
  wilcox.test(gene_expr_for_resp,gene_expr_for_nonresp)
pv_w[i] <-t2$p.value
}
}

pv
pv_w

# поправка на значения 
pv_adjusted <- p.adjust(pv, method="fdr")


library(ggplot2)
number_of_samples<-ncol(aggregated_expr)

df_long_for_2_genes<- data.frame("gene_name_column"=c(rep("A2M",number_of_samples),rep("A2LD1",number_of_samples)),
                            "expression"=c(unlist(aggregated_expr["A2M",]),unlist(aggregated_expr["A2LD1",])))

# 2 gene
boxplot(
  expression ~ gene_name_column,
  data = df_long_for_2_genes)

ggplot(data=df_long_for_2_genes, aes(x = expression,color=gene_name_column)) + geom_histogram()+theme_classic()
ggplot(data=df_long_for_2_genes, aes(x = expression,fill = gene_name_column)) + geom_histogram()+theme_classic()

ggplot(data=df_long_for_2_genes, aes(y = expression,fill = gene_name_column)) + geom_boxplot()+theme_classic()

ggplot(data=df_long_for_2_genes, aes(x = expression,fill = gene_name_column)) + geom_histogram()+theme_classic()+ggtitle("My histogram") + xlab("expression")+ ylab("count")

# 3 gene
number_of_samples<-ncol(aggregated_expr)
df_long_for_3_genes<- data.frame("gene_name_column"=c(rep("A2M",number_of_samples),rep("A2LD1",number_of_samples),rep("7A5",number_of_samples)),
"expression"=c(unlist(aggregated_expr["A2M",]),unlist(aggregated_expr["A2LD1",]),unlist(aggregated_expr["7A5",])))

library(ggplot2)
boxplot(
  expression ~ gene_name_column,
  data = df_long_for_3_genes)

ggplot(data=df_long_for_3_genes, aes(y = expression,fill = gene_name_column)) + geom_boxplot()+theme_classic()

ggplot(data=df_long_for_3_genes, aes(x = expression,fill = gene_name_column)) + geom_histogram()+theme_classic()+ggtitle("My histogram") + xlab("expression") 

ggplot(data=df_long_for_3_genes, aes(x=gene_name_column, y= expression, color="red"))+geom_boxplot()

ggplot(data=df_long_for_3_genes, aes(x=gene_name_column, y= expression, size = 0.25, alpha=0,1)) + geom_point(colour = "green") + ggtitle("Expressions genes") + xlab("gene_name_column")+ ylab("expression")
# aes(colour = factor(cyl),(shape = factor(cyl)

library(corrplot)               
M = cor(aggregated_expr,method = "s")
corrplot(M, method = 'color', order = 'alphabet')
corrplot(M, method = 'square', order = 'FPC', type = 'lower', diag = FALSE)
corrplot(M, method = 'square', diag = FALSE, order = 'hclust', addrect = 3, rect.col = 'white', rect.lwd = 3, tl.pos = 'd')


cor.test(unlist(aggregated_expr["7A5",]),unlist(aggregated_expr[ "A1BG",]),  method = "s")

# 
png("test_plot_Valeria.png")
corrplot.mixed(M, lower.col = "black", upper = "color", tl.col = "black", tl.cex = 0.8, tl.pos = "lt",
               order = 'hclust',
               hclust.method="ward.D2",
               diag = "u",number.cex = 0.5,
               sig.level = c(0.001, 0.01, 0.05), 
               pch.cex = 0.6,
               insig = 'label_sig',
               mar=c(0,1,0,0))
dev.off()


# extract metadata
md <- read.csv("C:/Users/rodic/OneDrive/?????????/GSE52219_series_matrix.txt",sep="\t",skip=28)[1:32,]
table(unlist(md[10,-1]))
number_of_patients<-ncol(md)-1
numbers_for_responders <- which(md[10,-1]=="response to mvac (downstaged to p0 or t1 w/ nac): yes")


numbers_for_nonresponders <- which(md[10,-1]=="response to mvac (downstaged to p0 or t1 w/ nac): no")
color_vector<-rep("",number_of_patients)
color_vector[numbers_for_responders]<-"green"
color_vector[numbers_for_nonresponders]<-"red"


color_vector

p<-prcomp(t(as.matrix(aggregated_expr[1:500,-1])))
pov<-p$sdev^2/sum(p$sdev^2)

library(scatterplot3d)
s3d <- scatterplot3d(p$x[,1], p$x[,2], p$x[,3], angle=60, pch = 7,type="p",xlab=paste0("PC1 (",round(pov[1]*100,digits = 2),"%)"),ylab=paste0("PC2 (",round(pov[2]*100,digits = 2),"%)"),zlab=paste0("PC3 (",round(pov[3]*100,digits = 2),"%)"),cex.symbols = 0.5, cex.axis = 0.7, cex.lab = 0.7, col = color_vector)


plot(p$x[,1], p$x[,2], pch = 18, type="p",xlab=paste0("PC1 (",round(pov[1]*100,digits = 2),"%)"),ylab=paste0("PC2 (",round(pov[2]*100,digits = 2),"%)"),cex.symbols = 0.5, cex.axis = 0.6, cex.lab = 0.8, col =color_vector, mar=c(1,2,1,1),main="2D PCA visualization")



# Heatmap
df_for_heatmap<-scale(as.matrix(aggregated_expr[1:20,]))
df = df_for_heatmap
heatmap(df, scale = "none")
color_pallete_rwb<- colorRampPalette(c("red", "white", "blue"))(256)

library("RColorBrewer")
heatmap(df, scale = "none", col = color_pallete_rwb,
        ColSideColors = color_vector)


library(ComplexHeatmap)
Heatmap(df, name = "Gene", #title of legend
        column_title = "Variables", row_title = "Samples",
        row_names_gp = gpar(fontsize = 8)) 


library("d3heatmap")
d3heatmap(scale(df), colors = "RdYlBu",
          k_row = 4,
          k_col = 2 )

d3heatmap(scale(df), colors = "RdYlBu")


  # ??????? ??????? 
library(tidyverse)
# aggregated_expr?
relig_income
relig_income %>%
  pivot_longer(!religion, names_to = "income", values_to = "count")

aggregated_expr$Gene_names <- rownames(aggregated_expr)


long_df_pivot <- aggregated_expr[1:3,]%>% #as_tibble%>%
  pivot_longer( !Gene_names, names_to = "Samples", values_to = "Expression")

# rm (df_) удалить из environment
# ctrl+f изменить название в нескольких местах


# Hc
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


