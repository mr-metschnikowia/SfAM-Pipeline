library(data.table)
library(plyr)
library(ggplot2)
library(RColorBrewer)
# importing prerequisites

df <- fread('C:/Users/Rhino/blast+/Outputs/output.txt',select=c(2,5))
# importing relevant columns from .csv table of blast output
krusk <- kruskal.test(df$qcovs,df$staxid)
chars <- capture.output(print(krusk))
writeLines(chars, con = file("C:/Users/Rhino/PycharmProjects/sfam/Figures/Krusk.txt"))
# data is not normally distributed and not matched pairs
# Kruskal-Wallis test is used to explore variation in qcovs between strains
# output is saved as .txt file
unstacked.qcovs = unstack(df[,c(2,1)])
png('C:/Users/Rhino/PycharmProjects/sfam/Figures/boxplot.png')
boxplot(unstacked.qcovs,col='yellow')
mtext('Strain Index', side=1, line=3, col='blue')
mtext('Query Coverage', side=2, line=3, col='blue', las=0)
title(main='Box and whisker plot of query coverage across yeast strains', col='black',cex=1.2)
dev.off()
# boxplot is created (in .png file) to visualise data

df2 <- fread('C:/Users/Rhino/PycharmProjects/sfam/pandas_dfs/df1.csv',select=c(3,4,9))
# importing relevant columns from .csv table of blast output
df3 <- ddply(df2,~staxid+clusters,summarise,accession_count=length(unique(sacc)))
df3$staxid <- as.factor(df3$staxid)
# creating new dataframe: number of unique accesion codes counted, aggregated on tax id and cluster
png('C:/Users/Rhino/PycharmProjects/sfam/Figures/barchart.png')
myColors <- brewer.pal(5, "Set3")
names(myColors) <- length(unique(df3$sacc))
bar <- ggplot(df3, aes(x=clusters, y=accession_count, fill=staxid)) + geom_bar(stat="identity") + scale_colour_manual(values=myColors)
bar <- bar + ggtitle("Distribution of gene clusters across accessions and across strains") + theme(plot.title=element_text(face="bold"))
bar <- bar + scale_fill_discrete(name="Strain\nIndex")
bar <- bar + xlab("Cluster") + ylab("Accession Count")
bar
dev.off()
# barchart is created (in .png file) to visualise data

