#install packages("devtools")
library(devtools)
#devtools::install_github('zdk123/SpiecEasi')
#install packages("SpiecEasi")
library(SpiecEasi)
#install packages("igraph")
library(igraph)

remove__ <- function(s, patten='__') {
  if (grepl(patten, s)) {
    s <- strsplit(s, patten)[[1]][2]
  }
  return(s)
}

setwd("E:/desktop/networkb-f/ck")
otutab <- read.table("ck.txt", sep = '\t', row.names = 1, header = T)
norm = t(t(otutab)/colSums(otutab,na=T))*100
idx = rowMeans(norm) > 0.00000000001
otutab_sub = otutab[idx,]
otu <- as.data.frame(t(otutab_sub))
set.seed(123)
otu.sparcc <- sparcc(otu, iter = 20, inner_iter = 10, th = 0.1)
sparcc0 <- otu.sparcc$Cor  
colnames(sparcc0) <- colnames(otu)
rownames(sparcc0) <- colnames(otu) 
write.table(sparcc0, 'sparcc0.txt', sep = '\t', col.names = NA, quote = FALSE)
n = 40
for (i in 1:n) {
  otu.boot <- sample(otu, replace = TRUE)  #bootstrap
  otu.sparcc_boot <- sparcc(otu.boot, iter = 20, inner_iter = 10, th = 0.1)  #sparcc 参数设置和上文保持一致
  sparcc_boot <- otu.sparcc_boot$Cor
  colnames(sparcc_boot) <- colnames(otu.boot)
  rownames(sparcc_boot) <- colnames(otu.boot)
  write.table(sparcc_boot, paste('sparcc_boot', i, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  
}

p <- sparcc0
p[p!=0] <- 0

for (i in 1:n) {
  p_boot <- read.delim(paste('sparcc_boot', i, '.txt', sep = ''), sep = '\t', row.names = 1)
  p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
}

p <- p/n
write.table(p, 'pvals.two_sided.txt', sep = '\t', col.names = NA, quote = FALSE)
corr <- as.data.frame(sparcc0)
pval <- as.data.frame(p)
corr <- corr[match(rownames(pval), rownames(corr)), match(rownames(pval), colnames(corr)), drop=F]
edges <- data.frame(from=rownames(corr)[row(corr)[upper.tri(corr)]], to=colnames(corr)[col(corr)[upper.tri(corr)]], corr=corr[upper.tri(corr)], pval=pval[upper.tri(pval)], stringsAsFactors=F)
edges <- edges[edges$pval < 0.01 & abs(edges$corr) > 0.8, ]

edges$type <- ifelse(edges$corr>0, 'positive', 'negative')
write.table(edges, paste0('network_edge.xls'), sep='\t', row.names=F, quote=F)

g <- graph.data.frame(d=edges, directed=FALSE)
# tname <- names(degree(g)[degree(g) > 0])
# if(length(tname) > 0) g <- induced_subgraph(g, tname)

cent <- data.frame(degree=degree(g), betweenness=betweenness(g), closeness= closeness(g), eigen=eigen_centrality(g)$vector, subgraph=subgraph_centrality(g))
write.table(cent, paste0('network_co-occurrence_node_infor.xls'), sep='\t', col.names=NA, quote=F)

sink(paste0('network_co-occurrence_network_property.xls'))
cat(sprintf('Number of edges\t%d\n',length(E(g))))
cat(sprintf('Number of vertices\t%d\n',length(V(g))))
cat(sprintf('Connectance\t%.5f\n',edge_density(g,loops=FALSE)))
cat(sprintf('Average degree\t%.5f\n', mean(degree(g))))
cat(sprintf('Average path length\t%.5f\n', average.path.length(g)))
cat(sprintf('Diameter\t%d\n', diameter(g, directed=FALSE, unconnected=TRUE, weights=NULL)))
cat(sprintf('Edge connectivity\t%d\n', edge_connectivity(g)))
cat(sprintf('Clustering coefficient\t%.5f\n', transitivity(g)))
cat(sprintf('Betweenness centralization\t%.5f\n', centralization.betweenness(g)$centralization))
cat(sprintf('Degree centralization\t%.5f\n', centralization.degree(g)$centralization))
sink(NULL)

degree <- degree(g)
nodes <- data.frame(name=names(degree), degree, label=sapply(names(degree), remove__))

if (!is.null("taxonomy.txt")){
  class <- read.table("taxonomy.txt", header=T, sep='\t', check.names=F, stringsAsFactors=F, comment.char='', row.names=1)
  class$otu <- rownames(class)
  
  color <- sapply(rownames(nodes), function(x) ifelse(x %in% class[,"otu"], strsplit(as.character(class[match(x, class[, "otu"]), "phylum"]), '__')[[1]][2], 'Other'))
  # color <- sapply(color, function(x) ifelse(length(which(color==x)) > 10, x, 'Other'))
  if ('Other' %in% color) {
    unother <- unique(color[color!='Other'])
    color <- factor(color, levels=c(unother, 'Other'))
  }
  nodes$color <- color
}
write.table(nodes, paste0('network_node.xls'), sep='\t', row.names=F, quote=F)

g <- graph_from_data_frame(edges, directed = FALSE, nodes)
write_graph(g, paste0('network_co-occurrence.graphml'), 'graphml')
