# adaptor: blastx is used to find alignment between translated pul genes and proteins > DNA sequences are obtained from protein sequences for each strain

import os
# importing prerequisites

def gene_dna_from_protein(species):
    # os.system("cd /home/centos/blast+")
    os.system("makeblastdb -in /home/centos/blast+/{}/protein.faa -out /home/centos/blast+/dbs/{}_proteins -parse_seqids -dbtype prot".format(species,species))
    # make database from protein faa
    # os.system("cd /home/centos/blast+")
    os.system('blastx -query /home/centos/project/genes/batch.txt -db {}_proteins -out /home/centos/blast+/outputs/{}_proteins.txt -outfmt "10 qacc stitle qcovs pident" -qcov_hsp_perc 70'.format(species,species))
    # blast pul gene dna against protein database of other strain
    locus_tags = []
    with open(r'/home/centos/blast+/outputs/{}_proteins.txt'.format(species), 'r') as f:
        data = f.readlines()
        for line in data:
            target = 'locus_tag='
            new_line = line[line.find(target):]
            new_line = new_line[len(target): new_line.find(']')]
            locus_tags.append(new_line)
    # locus_tags of hits are extracted
    with open(r'/home/centos/blast+/{}/genes.fna'.format(species), 'r') as f:
        data = f.read()
        data = data.split('>')
        # each gene is item in list
        genes = []
        for tag in locus_tags:
            for gene in data:
                result = gene.find(tag)
                if result > -1:
                    new_gene = '>' + gene
                    genes.append(new_gene)
    with open(r'/home/centos/project/genomes/{}_genes.txt'.format(species),'w') as f:
        f.writelines(genes)
    # dna sequences for each locus tag are found in .fna file of cds

if __name__ == '__main__':
    gene_dna_from_protein('k.lactis')
    gene_dna_from_protein('c.auris')
    gene_dna_from_protein('apc1.2')
