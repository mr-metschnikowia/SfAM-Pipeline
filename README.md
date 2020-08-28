# SfAM-Pipeline
Pipeline optimised to run on Test VM hosted on CLIMB.

User Input: 
1. .fasta files containing query sequences (/home/centos/project/genes)
2. .fasta files containing subject genomes (/home/centos/project/genomes)
3. Taxonomic map (accession-organismID)(/home/centos/blast+/tax_map/tax_map1.txt)
4. Predicted clusters (based on histogram)

IMPORTANT: At the moment, for some reason, file data is not overwritten but added to, therefore it is necessary to either delete all existing output files or change names of files before each run

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

All figures stored in - /home/centos/project/outputs
blast+ output - /home/centos/blast+/outputs
