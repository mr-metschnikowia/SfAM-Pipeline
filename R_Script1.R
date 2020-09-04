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

df2 <- fread('/home/centos/project/dfs/df1.csv',select=c(3,4,11))
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

df4 <- fread('/home/centos/project/dfs/df1.csv',select=c(2,3,12))
png('/home/centos/project/outputs/hist2.png')
hist(df4$`dnds`,xlab='dn/ds',main='Distribution of dn/ds')
# histogram of dn/ds plotted
dev.off()
df7 <- aggregate(df4$dnds, by=list(staxid=df4$staxid, qacc=df4$qacc), FUN=mean)
df7 <- rename(df7, c("x"="mean_dnds"))
png('/home/centos/project/outputs/barchart2.png')
no_colours = length(unique(df7$qacc))
myColors <- brewer.pal(no_colours, "Set1")
bar2 <- ggplot(df7, aes(x=staxid, y=mean_dnds, fill=qacc)) + geom_bar(stat="identity") + scale_colour_manual(values=myColors)
bar2 <- bar2 + ggtitle("Distrbituion of mean dn/ds across genes and across strains") + theme(plot.title=element_text(face="bold"))
bar2 <- bar2 + scale_fill_discrete(name="Gene")
bar2 <- bar2 + xlab("Strain") + ylab("mean dn/ds")
bar2
dev.off()
# multivariate bar chart of dn/ds by strain and genes

df5 <- fread('/home/centos/project/dfs/df1.csv',select=c(3,13))
df5 <- aggregate(df5$intergenic_space, by=list(staxid=df5$staxid), FUN=sum)
df5 <- rename(df5, c("x"="intergenic_space"))
df6 <- ddply(df2,~staxid,summarise,accession_count=length(unique(sacc)))
JoinedDT <- merge(df5,df6)
JoinedDT$staxid <- as.factor(JoinedDT$staxid)
png('/home/centos/project/outputs/hybrid.png')
scale <- max(JoinedDT$intergenic_space)/max(JoinedDT$accession_count)
hybrid <- ggplot(JoinedDT) + 
  geom_col(aes(x = staxid, y = intergenic_space), size = 1, color = "darkblue", fill = "white") +
  geom_line(aes(x = staxid, y = scale*accession_count), size = 1.5, color="red", group = 1) + 
  scale_y_continuous(sec.axis = sec_axis(~./scale, name = "Accession Count"))
hybrid <- hybrid + ggtitle("Distribution of intergenic space across strains and accessions") + theme(plot.title=element_text(face="bold"))
hybrid <- hybrid + xlab("Strain") + ylab("Total Intergenic Space")
hybrid
dev.off()
# hybrid (bar/line) of intergenic space by strain and number of accessions
