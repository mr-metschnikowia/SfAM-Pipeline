import os
import pandas as pd
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
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

def make_blast_db():
    os.system('makeblastdb -in /home/centos/project/genomes/custom.txt -dbtype nucl -parse_seqids -out /home/centos/blast+/dbs/custom_db -taxid_map /home/centos/blast+/tax_maps/tax_map1.txt')
    # interacts with command line - blastable database created from 'pre-database'

def batch_blastn():
    os.system("export BLASTDB='/home/centos/blast+/dbs'")
    os.system('blastn -query /home/centos/project/genes/batch.txt -db custom_db -out /home/centos/blast+/outputs/output.txt -outfmt "10 qacc staxid sacc pident qcovs sstart send"')
    # runs blast+ (blastn): batch query is blasted against custom database > output.txt file is pooped out
    with open(r'/home/centos/blast+/outputs/output.txt','r') as f:
        lines = f.readlines()
    lines.insert(0,'qacc,staxid,sacc,pident,qcovs,sstart,send\n')
    with open(r'/home/centos/blast+/outputs/output.txt', 'w') as f:
        f.writelines(lines)
    # title line is inserted

def dataframe():
    global df
    df = pd.read_csv(r'/home/centos/blast+/outputs/output.txt')
    df = df.sort_values(by=['qcovs'])
    # blast+ output is read into pandas dataframe
    # dataframe is ordered based on qcovs

def hist_me_up():
    df.hist(column = 'qcovs')
    plt.ylabel('frequency'), plt.xlabel('query coverage'), plt.title('Distribution of query coverage')
    plt.savefig(r'/home/centos/project/outputs/histogram.png')
    # histogram of qcovs is created
    # histogram is stored as figure1 in figures folder

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

def R_script():
    df.to_csv(r'/home/centos/project/dfs/df1.csv')
    # dataframe exported as .csv file
    os.system('/usr/bin/Rscript /home/centos/project/code/Script1.R')
    # R script 'Script1.R' is run

if __name__ == '__main__':
    pre_db()
    batch_q()
    make_blast_db()
    batch_blastn()
    dataframe()
    hist_me_up()
    clusterize()
    R_script()
    # calling functions
