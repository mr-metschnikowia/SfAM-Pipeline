# SfAM-Pipeline
User Input: 
1. .fasta files containing query sequences
2. .fasta files containing subject genomes
3. File path adjustments
4. Predicted clusters (based on histogram)

Prerequisites:
1. Relies on R script: https://github.com/mr-metschnikowia/Pipeline-Mod3

Automation:
1. .fasta files containing query sequences are concatenated to form batch query
2. Custom blastable database created from .fasta files containing subject genomes
3. Interaction with command line facilitates batch homology search with custom output (qacc,staxid,sacc,pident,qcovs,sstart,send)
4. Ouput is converted to pandas dataframe 
5. Histogram of qcovs is plotted to show distribution (ouput: .png file)
6. Hits are clustered based on qcovs using kmeans
7. Updated dataframe is exported as .csv file and relevant columns are read into R script
8. Kruskal-Wallis test conducted on qcovs aggregated by Strain (output: .txt file)
9. Boxplot of qcovs distribution for each strain plotted (output: .png file)
10. New metric: Number of unique accessions counted for each cluster aggregated on strain id
11. Multivariable bar chart created from data (output: .png file)
