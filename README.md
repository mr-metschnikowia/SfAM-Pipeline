# SfAM-Pipeline
Pipeline optimised to run on Test VM hosted on CLIMB.

User Input: 
1. .fasta files containing query sequences (/home/centos/project/genes)
2. .fasta files containing subject genomes (/home/centos/project/genomes)
3. Predicted clusters (based on histogram)
4. File path for databases - '/home/centos/blast+/dbs'

IMPORTANT: For clearance after each run use file_clearance.py

Automation:
1. taxonomic map is created from accessions of interest
2. .fasta files containing query sequences are concatenated to form batch query
3. Custom blastable database created from .fasta files containing subject genomes
4. Interaction with command line facilitates batch homology search with custom output (qacc,staxid,sacc,pident,qcovs,sstart,send)
5. Ouput is converted to pandas dataframe 
6. Histogram of qcovs is plotted to show distribution (ouput: .png file)
7. Hits are clustered based on qcovs using kmeans
8. master dataframe is fractured into gene specific dataframes
9. Gene specific dataframes are exported as .csv files
10. Gene specific dataframes are used for the following:
      Kruskal-Wallis is conducted to test for significant difference in query coverage across strains (seperate test for each gene)
      Boxplots are created reflecting distribution of query coverage for each strain (separate set for each gene) 
      For each gene: distribution of gene clusters across accessions in each strain 
11. Gene specific dataframes are further fractured into gene-cluster specific dataframes
12. DNDS is performed for each of these new dataframes - new column added for dnds of each hit against all other hits 
13 PLOT DNDS???????????????
14. master dataframe is split into new dataframes based on cluster
15. Intergenic sequences are examined separatly for each cluster:
      Length of intergenic sequences between genes on same accession in same strain is calculated
      Each length is assigned a standardised gene pair
16. Cluster dataframes with intergenic information are read into R
17. Multivariate bar chart is plotted to relfect distribution of intergenic space across gene pairs in different strains (cluster specific):
      For this a new metric is created: Average intergenic space per accession (to take accession count into account)
18. Specific intergenic sequence similarity???????

All figures stored in - /home/centos/project/outputs
blast+ output - /home/centos/blast+/outputs
