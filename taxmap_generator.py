# tax_map generator

import os

def make_blast_db1(input,output):
    os.system('makeblastdb -in /home/centos/project/genomes/{} -dbtype nucl -parse_seqids -out /home/centos/blast+/dbs/{}'.format(input,output))
    # interacts with command line - blastable database created from 'pre-database'

def make_blast_db2(input,output,tax_map):
    os.system('makeblastdb -in /home/centos/project/genomes/{} -dbtype nucl -parse_seqids -out /home/centos/blast+/dbs/{} -taxid_map /home/centos/blast+/tax_maps/{}'.format(input,output,tax_map))
    # interacts with command line - blastable database created from 'pre-database'

def batch_blastn(query,db,output):
    os.system("export BLASTDB='/home/centos/blast+/dbs'")
    os.system('blastn -query /home/centos/project/genes/{} -db {} -out /home/centos/blast+/pre_taxmaps/{} -outfmt "10 sacc"'.format(query,db,output))
    # runs blast+ (blastn): batch query is blasted against custom database > output.txt file is pooped out

if __name__ == '__main__':
    tax_dbs = []
    x = 0
    for file in os.listdir(r'/home/centos/project/genomes'):
        database = '{}'.format(x)
        make_blast_db1(file, database)
        tax_dbs.append(database)
        x += 1
    # database made from each genome (FASTA)
    x = 0
    tax_dbs = [0,1,2,3,4]
    for db in tax_dbs:
        y = '{}.txt'.format(x)
        batch_blastn('batch.txt', db, y)
        x += 1
    # batch query blasted against each genome to identify accessions of interest
    x = 0
    for file in os.listdir(r'/home/centos/blast+/pre_taxmaps'):
        path = r'/home/centos/blast+/pre_taxmaps/{}'.format(file)
        with open(path, 'r') as f:
            data = f.readlines()
        new_data = []
        for line in data:
            no_end = line.strip()
            new_line = no_end + ' {}\n'.format(x)
            new_data.append(new_line)
        x += 1
        with open(r'/home/centos/blast+/tax_maps/tax_map1.txt','a') as g:
            g.writelines(new_data)
    # taxonomic map is created with all accessions of interest in same strain assigned same taxonomic id
