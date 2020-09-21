# SfAM-Pipeline
Pipeline optimised to run on Test VM hosted on CLIMB.

User Input: 
1. .fasta files containing query sequences (/home/centos/project/genes)
2. .fasta files containing subject genomes (/home/centos/project/genomes) - names should be strain name (no '.'s in name)
3. Predicted clusters (based on histogram)
4. File path for databases - '/home/centos/blast+/dbs'

IMPORTANT: For clearance after each run use file_clearance.py

Automation:
1. get_dna.py: In each strain pul gene homologs are found using one of two methods - A protein homology search, followed by extraction of corresponding DNA sequences from .fna cds file, or (closely related species only) through a simple batch balstn homology search against genomes
2. needleall.py: Global alignments are performed using needleall EMBOSS between APC1.2 pul genes and pul genes from other species/strains.
3. Needleall output is processed to facilitate further analyses.
4. Pipeline.py: 
> historgram to show distribution of needle similarity score (normalised) is created
> kmeans clustering is used to assign each hit to a cluster based on similarity score
> hits are segregated on genes
> dnds is calculated for each hit versus all other hits, segregated by gene and cluster
> mean dnds is calculated for each hit
> intergenic space between each gene pair as well as total cluster intergenic space is calculated for each accession in each strain, segregated by cluster
5. R_Script1:
> boxplot showing distribution of similarity score for each gene is created
> histogram for each gene in each cluster is created to show distribution of dn/ds
> multivariate bar charts are created to show distribution of intergenic space between gene pairs across strains and variation in total cluster size between strains
> boxplot is also created to reflect distribution of intergenic space for each gene pair in each cluster
> relationship between gene length and mean dnds is plotted and a tested for correlation (Pearson's/Spearman's)

All figures stored in - /home/centos/project/outputs
blast+ output - /home/centos/blast+/outputs
