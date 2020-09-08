# SfAM-Pipeline
Pipeline optimised to run on Test VM hosted on CLIMB.

User Input: 
1. .fasta files containing query sequences (/home/centos/project/genes)
2. .fasta files containing subject genomes (/home/centos/project/genomes)
3. Predicted clusters (based on histogram)
4. File path for databases - '/home/centos/blast+/dbs'

IMPORTANT: At the moment, for some reason, file data is not overwritten but added to, therefore it is necessary to either delete all existing output files or change names of files before each run

Automation:
1. .fasta files containing query sequences are concatenated to form batch query
2. Custom blastable database created from .fasta files containing subject genomes
3. Interaction with command line facilitates batch homology search with custom output (qacc,staxid,sacc,pident,qcovs,sstart,send)
4. Ouput is converted to pandas dataframe 
5. Histogram of qcovs is plotted to show distribution (ouput: .png file)
6. Hits are clustered based on qcovs using kmeans
7. dn/ds is calculated for each hit relative to APC1.2
8. Updated dataframe is exported as .csv file and relevant columns are read into R script
9. Kruskal-Wallis test conducted on qcovs aggregated by Strain (output: .txt file)
10. Boxplot of qcovs distribution for each strain plotted (output: .png file)
11. New metric: Number of unique accessions counted for each cluster aggregated on strain id
12. Multivariable bar chart created from data (output: .png file)
13. Histogram of dn/ds is plotted
14. Multivariate barchart plotted reflecting mean dn/ds per gene per strain 
15. Number of chromosomes over which homolgs are distributed calculated for each strain
16. Total amount of intergenic space between homologs on same accession in each strain calculated
17. Hybrid(line/bar) created reflecting distribution of intergenic space between strains and accessions
18. All intergenic sequences are extracted
19. blastn search is conducted with intergenic sequences from APC1.2 as the query and intergenic sequences from other strains as subjects
20. Histogram of query coverage plotted for intergenic sequences
21. Bar chart plotted showing mean query coverage for intergenic sequences across strains
22. 'N's are accounted for in intergenic sequence comparison - graph plotted to reflect number of ambiguous sequences per strain 
23. taxonomic map is created from accessions of interest

All figures stored in - /home/centos/project/outputs
blast+ output - /home/centos/blast+/outputs
