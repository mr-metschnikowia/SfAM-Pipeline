# adaptor: blastx is used to find alignment between translated pul genes and proteins > DNA sequences are obtained from protein sequences for each strain

import os
import pandas as pd
# importing prerequisites

def gene_dna_from_genome(file_extension, species):
    os.system("makeblastdb -in /home/centos/blast+/{}/genome{} -out /home/centos/blast+/dbs/{}_genome -parse_seqids -dbtype nucl".format(species, file_extension, species))
    # make database from genome
    os.system('blastn -query /home/centos/project/genes/batch.txt -db {}_genome -out /home/centos/blast+/outputs/blastn_output/{}_genes.txt -outfmt "10 sacc sstart send sseq"'.format(species, species))
    # obtain dna sequences for pul homologs
    with open('/home/centos/blast+/outputs/blastn_output/{}_genes.txt'.format(species), 'r') as f:
        lines = f.readlines()
    lines.insert(0, 'sacc,sstart,send,sseq\n')
    with open('/home/centos/blast+/outputs/blastn_output/{}_genes.txt'.format(species), 'w') as f:
        f.writelines(lines)
    # title line is inserted

def process_blastn_data():
    for file in os.listdir('/home/centos/blast+/outputs/blastn_output/'):
        path = '/home/centos/blast+/outputs/blastn_output/' + file
        df = pd.read_csv(path)
        genes = []
        x = 0
        for i in range(len(df.sacc)):
            gene = '>' + df.sacc[i] + '.' + str(x) + ' [location={}..{}]'.format(df.sstart[i],df.send[i]) + ' [gbkey=CDS]\n' + df.sseq[i] + '\n'
            genes.append(gene)
            x += 1
        x = 0
        with open('/home/centos/project/genomes/{}'.format(file),'w') as f:
            f.writelines(genes)
    # blastn data converted into standardised form ready for emboss-needleall

def gene_dna_from_protein(species):
    os.system("makeblastdb -in /home/centos/blast+/{}/protein.faa -out /home/centos/blast+/dbs/{}_proteins -parse_seqids -dbtype prot".format(species,species))
    # make database from protein faa
    os.system('blastx -query /home/centos/project/genes/batch.txt -db {}_proteins -out /home/centos/blast+/outputs/{}_proteins.txt -outfmt "10 qacc stitle qcovs pident" -qcov_hsp_perc 60'.format(species,species))
    # blast pul gene dna against protein database of other strain
    with open(r'/home/centos/blast+/outputs/{}_proteins.txt'.format(species),'r') as f:
        lines = f.readlines()
    targets = []
    for line in lines:
        count = line.count(',')
        if float(line.split(',')[count].strip()) < 30:
            targets.append(lines.index(line))
    targets.sort(reverse=True)
    for target in targets:
        del lines[target]
    with open(r'/home/centos/blast+/outputs/{}_proteins.txt'.format(species),'w') as f:
        f.writelines(lines)
    # removing hits with pident < 30
    locus_tags = []
    with open(r'/home/centos/blast+/outputs/{}_proteins.txt'.format(species), 'r') as f:
        data = f.readlines()
        for line in data:
            if species.split('.')[0] == 's':
                target = 'protein_id='
            else:
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
    gene_dna_from_genome('.fsa_nt','ucd127')
    gene_dna_from_genome('.fasta','mp5')
    gene_dna_from_genome('.txt','bath1')
    gene_dna_from_genome('.fsa_nt','ap47')
    gene_dna_from_genome('.fna','277')
    process_blastn_data()
    gene_dna_from_protein('k.lactis')
    gene_dna_from_protein('c.auris')
    gene_dna_from_protein('apc1.2')
    gene_dna_from_protein('k.dobzhanskii')
    gene_dna_from_protein('z.mrakii')
    gene_dna_from_protein('k.pastoris')
    gene_dna_from_protein('l.thermotolerans')
    gene_dna_from_protein('k.marxianus')
    gene_dna_from_protein('s.eubayanus')
    gene_dna_from_protein('s.kudriavzevii')
    gene_dna_from_protein('s.paradoxus')
    gene_dna_from_protein('s.cerevisiae')
    gene_dna_from_protein('h.burtonii')
