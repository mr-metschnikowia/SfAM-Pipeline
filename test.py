# adaptor: blastx is used to find alignment between translated pul genes and proteins > DNA sequences are obtained from protein sequences
#          > DNA seqs added to database > pul gene BLAST against database
#          > hits integrated into master table?

import subprocess
# importing prerequisites

def make_blast_db1():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    proc.stdin.write(b"cd c:/users/rhino/blast+\n")
    proc.stdin.write(b"makeblastdb -in c:/users/rhino/blast+/k.lactis/protein.faa -out c:/users/rhino/blast+/blastdb/k.lactis_proteins -parse_seqids -dbtype prot\n")
    proc.stdin.close()
    proc.wait()
    # make database from protein faa

def blastx():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    proc.stdin.write(b"cd c:/users/rhino/blast+\n")
    proc.stdin.write(b'blastx -query c:/users/rhino/blast+/q_stage2/batch.txt -db k.lactis_proteins -out c:/users/rhino/blast+/Outputs/k.lactis_proteins.txt -outfmt "10 qacc stitle qcovs"\n')
    proc.stdin.close()
    proc.wait()
    # blast pul gene dna against protein database of other strain

def parse_protein_data():
    global locus_tags
    with open(r'c:/users/rhino/blast+/Outputs/k.lactis_proteins.txt', 'r') as f:
        data = f.readlines()
        locus_tags = []
        for line in data:
            target = 'locus_tag='
            new_line = line[line.find(target) + len(target): line.find(']')]
            locus_tags.append(new_line)
    # locus_tags of hits are extracted

def find_dna_sequences():
    with open(r'C:\Users\Rhino\blast+\k.lactis\genes.fna', 'r') as f:
        data = f.read()
        data = data.split('>')
        # each gene is item in list
        genes = []
        for tag in locus_tags:
            for gene in data:
                result = gene.find(tag)
                if result > -1:
                    gene = '>' + gene
                    genes.append(gene)
    with open(r'c:\users\rhino\blast+\db_stage1\k.lactis_genes.txt','w') as f:
        f.writelines(genes)
    # dna sequences for each locus tag are found in .fna file of cds

def make_blast_db2():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    proc.stdin.write(b"cd c:/users/rhino/blast+\n")
    proc.stdin.write(b"makeblastdb -in c:/users/rhino/blast+/db_stage1/k.lactis_genes.txt -out c:/users/rhino/blast+/blastdb/k.lactis_genes -parse_seqids -dbtype nucl\n")
    proc.stdin.close()
    proc.wait()
    # database made from dna sequences (from protein hits) of other strain

def blastn():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    proc.stdin.write(b"cd c:/users/rhino/blast+\n")
    proc.stdin.write(b'blastn -query c:\users\\rhino\\blast+\q_stage2\\batch.txt -db k.lactis_genes -out c:/users/rhino/blast+/Outputs/k.lactis_genes.txt -outfmt "10 qacc stitle qcovs"\n')
    proc.stdin.close()
    proc.wait()
    # pul gene dna blasted against database of gene dna (other strain)

if __name__ == '__main__':
    make_blast_db1()
    blastx()
    parse_protein_data()
    find_dna_sequences()
    make_blast_db2()
    blastn()
