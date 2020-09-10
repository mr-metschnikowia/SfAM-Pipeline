library(data.table)
library(plyr)
library(ggplot2)
library(RColorBrewer)
# importing prerequisites

reference <- fread ('/home/centos/project/dfs/reference.csv',select=c(11))
# creating reference

query_coverage_across_strains <- function(path,gene) {
  df <- fread(path,select=c(3,6))
  # importing relevant columns from .csv table of blast output
  krusk <- kruskal.test(df$qcovs,df$staxid)
  chars <- capture.output(print(krusk))
  path2 = sprintf("/home/centos/project/outputs/%s_krusk.txt",gene)
  writeLines(chars, con=file(path2))
  # data is not normally distributed and not matched pairs
  # Kruskal-Wallis test is used to explore variation in qcovs between strains
  # output is saved as .txt file
  unstacked.qcovs = unstack(df[,c(2,1)])
  boxplot(unstacked.qcovs,col='yellow')
  mtext('Strain Index', side=1, line=3, col='blue')
  mtext('Query Coverage', side=2, line=3, col='blue', las=0)
  title(main=gene, col='black',cex=1.2)
  boxplot
  # boxplot is created (in .png file) to visualise data
}

cluster_distribution <- function(path,gene) {
  df2 <- fread(path,select=c(3,4,11))
  # importing relevant columns from .csv table of blast output
  df3 <- ddply(df2,~staxid+clusters,summarise,accession_count=length(unique(sacc)))
  df3$staxid <- as.factor(df3$staxid)
  df3$clusters <- as.factor(df3$clusters)
  # creating new dataframe: number of unique accesion codes counted, aggregated on tax id and cluster
  path2 = sprintf('/home/centos/project/outputs/%s_cluster_distribution.png',gene)
  no_colours = length(unique(df3$clusters))
  myColors <- brewer.pal(no_colours, "Set3")
  bar <- ggplot(df3, aes(x=staxid, y=accession_count, fill=clusters)) + geom_bar(stat="identity") + scale_colour_manual(values=myColors)
  bar <- bar + ggtitle("Distribution of gene clusters across strains") + theme(plot.title=element_text(face="bold"))
  bar <- bar + scale_fill_discrete(name="Cluster")
  bar <- bar + xlab("Strain") + ylab("Accession Count")
  bar
  ggsave(path2)
  # barchart is created (in .png file) to visualise data
}

Intergenic_Space_Distribution <- function(cluster) {
  path = sprintf('/home/centos/project/dfs/cluster%s.csv',cluster)
  print(path)
  df5 <- fread(path,select=c(3,4,7,8))
  df5$staxid <- as.factor(df5$staxid)
  df6 <- ddply(df5,~gene_pair+staxid,summarise,accession_count=length(unique(sacc)))
  df7 <- aggregate(df5$intergenic_space, by=list(staxid=df5$staxid, gene_pair=df5$gene_pair), FUN=sum)
  df7 <- rename(df7, c("x"="intergenic_space"))
  JoinedDT <- merge(df6,df7)
  JoinedDT$average_intergenic_space_per_accession <- JoinedDT$intergenic_space/JoinedDT$accession_count
  path2 = sprintf('/home/centos/project/outputs/cluster%s_intergenic_space_distribution.png',cluster)
  no_colours = length(unique(JoinedDT$gene_pair))
  myColors <- brewer.pal(no_colours, "Set3")
  bar <- ggplot(JoinedDT, aes(x=staxid, y=average_intergenic_space_per_accession, fill=gene_pair)) + geom_bar(stat="identity") + scale_colour_manual(values=myColors)
  bar <- bar + ggtitle("Distribution of intergenic space across pul cluster in different strains") + theme(plot.title=element_text(face="bold"))
  bar <- bar + scale_fill_discrete(name="Gene\nPair")
  bar <- bar + xlab("Strain") + ylab("Total Intergenic Space/Accession")
  bar
  ggsave(path2)
  # multivariate bar chart plotted reflecting distribution of intergenic space across gene pairs across strains in cluster
}

png('/home/centos/project/outputs/boxplot.png')
par(mfrow = c(2,2), oma = c(0,0,2,0))
query_coverage_across_strains('/home/centos/project/dfs/pul1.csv','pul1')
query_coverage_across_strains('/home/centos/project/dfs/pul2.csv','pul2')
query_coverage_across_strains('/home/centos/project/dfs/pul3.csv','pul3')
query_coverage_across_strains('/home/centos/project/dfs/pul4.csv','pul4')
mtext("Graph to show query coverage (relative to APC1.2 - Strain 1) of homologs\n of each gene (pul1-4) in different yeast strains",side=3,outer=TRUE,padj=3, line=5)
dev.off()
par(mfrow = c(1,1))
cluster_distribution('/home/centos/project/dfs/pul1.csv','pul1')
cluster_distribution('/home/centos/project/dfs/pul2.csv','pul2')
cluster_distribution('/home/centos/project/dfs/pul3.csv','pul3')
cluster_distribution('/home/centos/project/dfs/pul4.csv','pul4')
for (val in unique(reference$clusters)) {
  try(Intergenic_Space_Distribution(val), silent=T)
}
