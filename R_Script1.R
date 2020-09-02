library(data.table)
library(plyr)
library(ggplot2)
library(RColorBrewer)
# importing prerequisites

df <- fread('/home/centos/blast+/outputs/output.txt',select=c(2,5))
# importing relevant columns from .csv table of blast output
krusk <- kruskal.test(df$qcovs,df$staxid)
chars <- capture.output(print(krusk))
writeLines(chars, con = file("/home/centos/project/outputs/Krusk.txt"))
# data is not normally distributed and not matched pairs
# Kruskal-Wallis test is used to explore variation in qcovs between strains
# output is saved as .txt file
unstacked.qcovs = unstack(df[,c(2,1)])
png('/home/centos/project/outputs/boxplot.png')
boxplot(unstacked.qcovs,col='yellow')
mtext('Strain Index', side=1, line=3, col='blue')
mtext('Query Coverage', side=2, line=3, col='blue', las=0)
title(main='Box and whisker plot of query coverage across yeast strains', col='black',cex=1.2)
dev.off()
# boxplot is created (in .png file) to visualise data

df2 <- fread('/home/centos/project/dfs/df1.csv',select=c(3,4,12))
# importing relevant columns from .csv table of blast output
df3 <- ddply(df2,~staxid+clusters,summarise,accession_count=length(unique(sacc)))
df3$staxid <- as.factor(df3$staxid)
# creating new dataframe: number of unique accesion codes counted, aggregated on tax id and cluster
png('/home/centos/project/outputs/barchart.png')
no_colours = length(unique(df3$clusters))
myColors <- brewer.pal(no_colours, "Set3")
names(myColors) <- length(unique(df3$sacc))
bar <- ggplot(df3, aes(x=clusters, y=accession_count, fill=staxid)) + geom_bar(stat="identity") + scale_colour_manual(values=myColors)
bar <- bar + ggtitle("Distribution of gene clusters across accessions and across strains") + theme(plot.title=element_text(face="bold"))
bar <- bar + scale_fill_discrete(name="Strain\nIndex")
bar <- bar + xlab("Cluster") + ylab("Accession Count")
bar
dev.off()
# barchart is created (in .png file) to visualise data

df4 <- fread('/home/centos/project/dfs/df1.csv',select=c(2,3,11))
png('/home/centos/project/outputs/hist2.png')
hist(df4$`dnds`,xlab='dn/ds',main='Distribution of dn/ds')
dev.off()
png('/home/centos/project/outputs/barchart2.png')
no_colours = length(unique(df4$qacc))
myColors <- brewer.pal(no_colours, "Set3")
bar2 <- ggplot(df4, aes(x=staxid, y=dnds, fill=qacc)) + geom_bar(stat="identity") + scale_colour_manual(values=myColors)
bar2 <- bar2 + ggtitle("Distribution of dn/ds across strains and across homologs") + theme(plot.title=element_text(face="bold"))
bar2 <- bar2 + scale_fill_discrete(name="Gene")
bar2 <- bar2 + xlab("Strain") + ylab("dn/ds")
bar2
dev.off()
