import os
import subprocess
import pandas as pd
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
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
    proc.stdin.write(b'makeblastdb -in C:\Users\Rhino\\blast+\db_stage2\custom.txt -dbtype nucl -parse_seqids -out C:\Users\Rhino\\blast+\\blastdb\custom_db -taxid_map C:\Users\Rhino\\blast+\\tax_map\\tax_map1.txt\n')
    proc.stdin.close()
    proc.wait()
    # Local BLASTable database is created from pre database

def batch_blastn():
    proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
    # facilitates command line interaction
    proc.stdin.write(b'blastn -query C:\Users\Rhino\\blast+\q_stage2\\batch.txt -db custom_db -out C:\Users\Rhino\\blast+\Outputs\output.txt -outfmt "10 qacc staxid sacc pident qcovs sstart send"\n')
    proc.stdin.close()
    proc.wait()
    # runs blast+ (blastn): batch query is blasted against custom database > output.txt file is pooped out
    with open(r'C:\Users\Rhino\\blast+\Outputs\output.txt','r') as f:
        lines = f.readlines()
    lines.insert(0,'qacc,staxid,sacc,pident,qcovs,sstart,send\n')
    with open(r'C:\Users\Rhino\\blast+\Outputs\output.txt', 'w') as f:
        f.writelines(lines)
    # title line is inserted

def dataframe():
    global df
    df = pd.read_csv(r'C:\Users\Rhino\blast+\Outputs\output.txt')
    df = df.sort_values(by=['qcovs'])
    # blast+ output is read into pandas dataframe
    # dataframe is ordered based on qcovs

def hist_me_up():
    df.hist(column = 'qcovs')
    plt.ylabel('frequency'), plt.xlabel('query coverage'), plt.title('Distribution of query coverage')
    plt.savefig(r'C:\Users\Rhino\PycharmProjects\sfam\Figures\\histogram.png')
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
    df.to_csv(r'C:\Users\Rhino\PycharmProjects\sfam\pandas_dfs\df1.csv')
    # dataframe exported as .csv file
    rscript = subprocess.call(['C:/Users/Rhino/Documents/R/R-3.6.2/bin/Rscript.exe', '--vanilla', 'C:/Users/Rhino/PycharmProjects/sfam/R_Scripts/Script1.R'], shell=True)
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
