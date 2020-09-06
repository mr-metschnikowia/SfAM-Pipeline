library(data.table)
library(plyr)
library(ggplot2)
# importing prerequisites

df <- fread('/home/centos/blast+/outputs/output2.txt',select=c(2,5))
#reading in relevant columns
df <- aggregate(df$qcovs, by=list(staxid=df$staxid), FUN=mean)
# taking mean of qcovs aggregated on staxid
df <- rename(df, c("x"="mean_qcovs"))
df$staxid <- as.factor(df$staxid)
png('/home/centos/project/outputs/bar3.png')
bar <- ggplot(df) + 
  geom_col(aes(x = staxid, y = mean_qcovs), size = 1, color = "red", fill = "white")
bar <- bar + ggtitle("Mean query coverage across yeast strains") + theme(plot.title=element_text(face="bold"))
bar <- bar + xlab("Strain") + ylab("Mean Query Coverage")
bar
dev.off()
# bar chart plotted of mean qcovs across strains and .png file is pooped out
