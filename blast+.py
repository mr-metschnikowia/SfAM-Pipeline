import os
import subprocess
# importing prereqs

def pre_db():
    master_string_1 = ''
    for file in os.listdir(r'C:\Users\Rhino\blast+\db_stage1'):
        with open(r'C:\Users\Rhino\blast+\db_stage1' + '\\' + file,'r') as f:
            contents = f.read()
            master_string_1 = master_string_1 + contents
    with open(r'C:\Users\Rhino\blast+\db_stage2\custom.txt', 'w') as f:
        f.write(master_string_1)
    # creates custom pre database from FASTA files in db_stage1 folder

def batch_q():
    master_string_2 = ''
    for file in os.listdir(r'C:\Users\Rhino\blast+\q_stage1'):
        with open(r'C:\Users\Rhino\blast+\q_stage1' + '\\' + file,'r') as f:
            contents = f.read() + '\n'
            master_string_2 = master_string_2 + contents
    with open(r'C:\Users\Rhino\blast+\q_stage2\batch.txt', 'w') as f:
        f.write(master_string_2)
    # creates batch query file from gene FASTA files

def make_blast_db():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    # facilitates command line interaction
    proc.stdin.write(b'cd C:\Users\Rhino\\blast+\n')
    proc.stdin.write(b'makeblastdb -in C:\Users\Rhino\\blast+\db_stage2\custom.txt -dbtype nucl -parse_seqids -out C:\Users\Rhino\\blast+\\blastdb\custom_db\n')
    proc.stdin.close()
    proc.wait()
    # Local BLASTable database is created from pre database

def batch_blastn():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    # facilitates command line interaction
    proc.stdin.write(b'blastn -query C:\Users\Rhino\\blast+\q_stage2\\batch.txt -db custom_db -out C:\Users\Rhino\\blast+\Outputs\output.txt -outfmt "10 qacc sacc pident qcovs sstart send"\n')
    proc.stdin.close()
    proc.wait()
    # runs blast+ (blastn): batch query is blasted against custom database > output.txt file is pooped out
    with open(r'C:\Users\Rhino\\blast+\Outputs\output.txt','r') as f:
        lines = f.readlines()
    lines.insert(0,'qacc,sacc,pident,qcovs,sstart,send\n')
    with open(r'C:\Users\Rhino\\blast+\Outputs\output.txt', 'w') as f:
        f.writelines(lines)
    # title line is inserted

if __name__ == '__main__':
    pre_db()
    batch_q()
    make_blast_db()
    batch_blastn()
    # calling functions
