import os
import pandas as pd
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
from dnds import dnds
# importing prereqs

def pre_db():
    master_string_1 = ''
    for file in os.listdir(r'/home/centos/project/genomes'):
        with open(r'/home/centos/project/genomes' + '/' + file,'r') as f:
            contents = f.read()
            master_string_1 = master_string_1 + contents
    with open(r'/home/centos/project/genomes/custom.txt', 'w') as f:
        f.write(master_string_1)
    # creates custom pre database from FASTA files in db_stage1 folder

def batch_q():
    master_string_2 = ''
    for file in os.listdir(r'/home/centos/project/genes'):
        with open(r'/home/centos/project/genes' + '/' + file,'r') as f:
            contents = f.read() + '\n'
            master_string_2 = master_string_2 + contents
    with open(r'/home/centos/project/genes/batch.txt', 'w') as f:
        f.write(master_string_2)
    # creates batch query file from gene FASTA files

def make_blast_db(input,output,tax_map):
    os.system('makeblastdb -in /home/centos/project/genomes/{} -dbtype nucl -parse_seqids -out /home/centos/blast+/dbs/{} -taxid_map /home/centos/blast+/tax_maps/{}'.format(input,output,tax_map))
    # interacts with command line - blastable database created from 'pre-database'

def batch_blastn(query,db,output):
    os.system("export BLASTDB='/home/centos/blast+/dbs'")
    os.system('blastn -query /home/centos/project/genes/{} -db {} -out /home/centos/blast+/outputs/{} -outfmt "10 qacc staxid sacc pident qcovs sstart send qseq sseq"'.format(query,db,output))
    # runs blast+ (blastn): batch query is blasted against custom database > output.txt file is pooped out
    with open(r'/home/centos/blast+/outputs/{}'.format(output),'r') as f:
        lines = f.readlines()
    lines.insert(0,'qacc,staxid,sacc,pident,qcovs,sstart,send,qseq,sseq\n')
    with open(r'/home/centos/blast+/outputs/{}'.format(output), 'w') as f:
        f.writelines(lines)
    # title line is inserted

def dataframe(path):
    global df
    df = pd.read_csv(path)
    # blast+ output is read into pandas dataframe

def hist_me_up(output):
    global empty_df
    empty_df = False
    try:
        df.hist(column='qcovs')
        plt.ylabel('frequency'), plt.xlabel('query coverage'), plt.title('Distribution of query coverage')
        plt.savefig(output)
    except ValueError:
        empty_df = True
        print('There is no data in this frame de data')
    # histogram of qcovs is created
    # histogram in figures folder
    # empty data frame is accomodated for

def clusterize():
    column = df.qcovs
    array = column.values
    prediction = int(input('Please predict number of clusters:'))
    compatible = array.reshape(-1,1)
    km = KMeans(n_clusters=prediction,random_state=0).fit(compatible)
    labels = km.labels_
    labels = [str(i) for i in labels]
    df['clusters'] = labels
    # 1D array created from dataframe column of interest
    # 2D array created for compatibility
    # k-means cluster analysis is performed on 2D array
    # cluster data is added to dataframe

def DNDS():
    column = []
    for i in range(len(df.qseq)):
        codons = ['-', 'K', 'M', 'B', 'V', 'S', 'W', 'D', 'Y', 'R', 'H']
        for codon in codons:
            while df.qseq[i].find(codon) > -1:
                df.sseq[i] = df.sseq[i][0: df.qseq[i].find(codon)] + df.sseq[i][df.qseq[i].find(codon) + 1:]
                df.qseq[i] = df.qseq[i][0: df.qseq[i].find(codon)] + df.qseq[i][df.qseq[i].find(codon) + 1:]
            while df.sseq[i].find(codon) > -1:
                df.qseq[i] = df.qseq[i][0: df.sseq[i].find(codon)] + df.qseq[i][df.sseq[i].find(codon) + 1:]
                df.sseq[i] = df.sseq[i][0: df.sseq[i].find(codon)] + df.sseq[i][df.sseq[i].find(codon) + 1:]
        # removing gaps
        while len(df.qseq[i]) % 3 > 0:
            df.qseq[i] = df.qseq[i][:-1]
            df.sseq[i] = df.sseq[i][:-1]
        # trimming sequences that are not multiples of 3
        seq_1 = df.qseq[i]
        seq_2 = df.sseq[i]
        try:
            x = round(dnds(seq_1, seq_2), 3)
        except ZeroDivisionError:
            x = 0
        column.append(x)
        # dnds calculated for each hit relative to APC1.2
    df['dnds'] = column
    # dnds column added to dataframe
    # x = 0
    # # y = 3
    # # seq = 'AAABBBCCCDD-'
    # # for i in range(int(len(seq)/3)):
    # #     if seq[x:y].find('-') > -1:
    # #         seq = seq[0:x] + seq[y:]
    # #     x = y
    # #     y = x + 3
    # # # alternative method for omitting unwanted codons, taking triplets into account

def intergenic_space_count():
    global df
    df.sort_values(['staxid', 'sacc', 'sstart'], inplace=True)
    df = df.reset_index(drop=True)
    # data ordered by staxid, sacc and sstart
    total_intergenic_space = []
    for i in range(len(df.staxid)):
        if df.sstart[i] > df.send[i]:
            x = df.sstart[i]
            y = df.send[i]
            df.sstart[i] = y
            df.send[i] = x
    # sstart and ssend are flipped if ssend < sstart for standardisation purposes
    for i in range(len(df.staxid) - 1):
        staxid = df.staxid[i]
        sacc = df.sacc[i]
        if staxid == df.staxid[i + 1]:
            if sacc == df.sacc[i + 1]:
                intergenic_space = df.sstart[i + 1] - df.send[i]
                total_intergenic_space.append(intergenic_space)
            else:
                intergenic_space = 0
                total_intergenic_space.append(intergenic_space)
        else:
            intergenic_space = 0
            total_intergenic_space.append(intergenic_space)
    total_intergenic_space.append(0)
    df['intergenic_space'] = total_intergenic_space
    # intergenic space between each gene on same chromosme/contig in same strain is calculated

def extract_intergenic_seqs():
    with open(r'/home/centos/project/genomes/custom.txt', 'r') as f:
        data = f.read()
    column = []
    for i in range(len(df.sacc)):
        if df.intergenic_space[i] > 0:
            data.find(df.sacc[i])
            chunk1 = data[data.find(df.sacc[i]):]
            chunk2 = chunk1[chunk1.find('\n'):]
            start = df.send[i] + 1
            end = df.sstart[i + 1]
            inter_seq = chunk2[start:end]
            column.append(inter_seq)
        else:
            inter_seq = 'N/A'
            column.append(inter_seq)
    df['inter_seqs'] = column
    # intergenic sequences are extracted and added to dataframe

def account_for_Ns():
    column = []
    for i in range(len(df.inter_seqs)):
        count = 0
        sequence = df.inter_seqs[i]
        if sequence == 'N/A':
            column.append(0)
        elif sequence.find('N') > -1:
            for codon in sequence:
                if codon == 'N':
                    count += 1
            df.intergenic_space[i] = df.intergenic_space[i] - count
            column.append(1)
        else:
            column.append(0)
    df['ambiguous_seqs'] = column
    # 'N's are accounted for: Removed from sequence length and number of ambigous sequences for is recorded for each strain

def batch_from_string():
    with open(r'/home/centos/project/genes/batch2.txt', 'w') as f:
        count = 1
        tax_codes = []
        for i in range(len(df.sacc)):
            if df.sacc[i] == 'CP034458' and type(df.inter_seqs[i]) == str:
                new_acc = '{}.{}'.format(df.sacc[i],count)
                tax_codes.append(new_acc + ' 1\n')
                f.write('>{}\n'.format(new_acc))
                f.write(df.inter_seqs[i])
                f.write('\n')
                count +=1
        with open(r'/home/centos/blast+/tax_maps/tax_map2.txt','w') as f:
            f.writelines(tax_codes)

def database_from_string():
    with open(r'/home/centos/project/genomes/custom2.txt', 'w') as f:
        scanned = []
        count = 1
        staxid = 1
        tax_codes = []
        for i in range(len(df.sacc)):
            if df.sacc[i] == 'CP034458':
                continue
            elif type(df.inter_seqs[i]) == str:
                if df.sacc[i] in scanned:
                    count += 1
                else:
                    count = 1
                    staxid += 1
                new_acc = '{}.{}'.format(df.sacc[i], count)
                tax_codes.append(new_acc + ' {}\n'.format(staxid))
                f.write('>{}\n'.format(new_acc))
                f.write(df.inter_seqs[i])
                f.write('\n')
                scanned.append(df.sacc[i])
        with open(r'/home/centos/blast+/tax_maps/tax_map2.txt','a') as f:
            f.writelines(tax_codes)

def R_script(script):
    os.system('/usr/bin/Rscript' + ' {}'.format(script))
    # R script is run

if __name__ == '__main__':
    batch_q()
    os.system('python3 /home/centos/project/code/taxmap_generator.py')
    # generating taxonomic map
    pre_db()
    make_blast_db('custom.txt','custom_db','tax_map1.txt')
    batch_blastn('batch.txt','custom_db','output.txt')
    dataframe(r'/home/centos/blast+/outputs/output.txt')
    hist_me_up(r'/home/centos/project/outputs/histogram.png')
    clusterize()
    DNDS()
    intergenic_space_count()
    extract_intergenic_seqs()
    account_for_Ns()
    df.to_csv(r'/home/centos/project/dfs/df1.csv')
    R_script('/home/centos/project/code/Script1.R')
    batch_from_string()
    database_from_string()
    make_blast_db('custom2.txt', 'custom2', 'tax_map2.txt')
    batch_blastn('batch2.txt', 'custom2', 'output2.txt')
    dataframe(r'/home/centos/blast+/outputs/output2.txt')
    hist_me_up(r'/home/centos/project/outputs/histogram3.png')
    if empty_df == False:
        R_script('/home/centos/project/code/Script2.R')
    # calling functions
    # data frames are exported as .csv files
    # Second R script is run if data frame contains data
