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
    global prediction
    column = df.qcovs
    array = column.values
    prediction = int(input('Please predict number of clusters:'))
    compatible = array.reshape(-1,1)
    km = KMeans(n_clusters=prediction,random_state=0).fit(compatible)
    labels = km.labels_
    labels = [str(i) for i in labels]
    df['clusters'] = labels
    df.to_csv(r'/home/centos/project/dfs/reference.csv')
    # 1D array created from dataframe column of interest
    # 2D array created for compatibility
    # k-means cluster analysis is performed on 2D array
    # cluster data is added to dataframe
    # dataframe is shipped as .csv to provide reference for R script

def DNDS():
    global target_frame2
    for j in range(len(target_frame2.sseq)):
        column = []
        for i in range(len(target_frame2.sseq)):
            reference = target_frame2.sseq[j]
            reference_length = len(reference)
            seq2 = target_frame2.sseq[i]
            seq2_length = len(seq2)
            if seq2_length > reference_length:
                seq2 = seq2[:reference_length]
            else:
                reference = reference[:seq2_length]
            # ensuring equal length
            for c in range(reference_length):
                codons = ['-', 'K', 'M', 'B', 'V', 'S', 'W', 'D', 'Y', 'R', 'H']
                for codon in codons:
                    while reference.find(codon) > -1:
                        seq2 = seq2[0: reference.find(codon)] + seq2[reference.find(codon) + 1:]
                        reference = reference[0: reference.find(codon)] + reference[reference.find(codon) + 1:]
                    while seq2.find(codon) > -1:
                        reference = reference[0: seq2.find(codon)] + reference[seq2.find(codon) + 1:]
                        seq2 = seq2[0: seq2.find(codon)] + seq2[seq2.find(codon) + 1:]
                # removing ambiguous codons
            while len(reference) % 3 > 0:
                reference = reference[:-1]
                seq2 = seq2[:-1]
            # trimming sequences that are not multiples of 3
            try:
                x = round(dnds(reference, seq2), 3)
            except ZeroDivisionError:
                x = 'N/A'
            except ValueError:
                x = 'math_error'
            column.append(x)
            # dnds calculated for each hit relative to all other hits
        column_name = str(target_frame2.staxid[j]) + '-{}'.format(target_frame2.sacc[j])
        target_frame2[column_name] = column
        # dnds column for each hit added to dataframe

def intergenic_space_count():
    global ready_df
    ready_df.sort_values(['staxid', 'sacc', 'sstart'], inplace=True)
    ready_df = ready_df.reset_index(drop=True)
    # data ordered by staxid, sacc and sstart
    total_intergenic_space = []
    for i in range(len(ready_df.staxid)):
        if ready_df.sstart[i] > ready_df.send[i]:
            x = ready_df.sstart[i]
            y = ready_df.send[i]
            ready_df.sstart[i] = y
            ready_df.send[i] = x
    # sstart and ssend are flipped if ssend < sstart for standardisation purposes
    for i in range(len(ready_df.staxid) - 1):
        staxid = ready_df.staxid[i]
        sacc = ready_df.sacc[i]
        if staxid == ready_df.staxid[i + 1]:
            if sacc == ready_df.sacc[i + 1]:
                intergenic_space = ready_df.sstart[i + 1] - ready_df.send[i]
                total_intergenic_space.append(intergenic_space)
            else:
                intergenic_space = 0
                total_intergenic_space.append(intergenic_space)
        else:
            intergenic_space = 0
            total_intergenic_space.append(intergenic_space)
    total_intergenic_space.append(0)
    ready_df['intergenic_space'] = total_intergenic_space
    # intergenic space between each gene on same chromosme/contig in same strain is calculated

def GenePair():
    global ready_df
    possible_pairs = ['pul1-pul1','pul1-pul2','pul1-pul3','pul1-pul4','pul2-pul2','pul2-pul3','pul2-pul4','pul3-pul3','pul3-pul4','pul4-pul4']
    column = []
    for i in range(len(ready_df.intergenic_space)):
        if ready_df.intergenic_space[i] > 0:
            gene_pair = '{}-{}'.format(ready_df.qacc[i],ready_df.qacc[i+1])
            if gene_pair not in possible_pairs:
                gene_pair = '{}-{}'.format(ready_df.qacc[i+1],ready_df.qacc[i])
                column.append(gene_pair)
            else:
                column.append(gene_pair)
        else:
            column.append('NA')
    ready_df['gene_pair'] = column
    # for each intergenic sequence standardised gene pair is allocated

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

def total_space():
    global ready_df
    global name
    fragments_of_total_area = []
    for i in range(len(ready_df.sstart)):
        gene_and_intergenic_space = (ready_df.send[i] - ready_df.sstart[i]) + ready_df.intergenic_space[i]
        fragments_of_total_area.append(gene_and_intergenic_space)
    ready_df['gene_length_and_intergenic_space'] = fragments_of_total_area
    # adds intergenic space to gene length
    staxids = []
    saccs = []
    total_space = []
    to_sum = []
    for i in range(len(ready_df.sstart)):
        if i < len(ready_df.sstart) - 1 and ready_df.sacc[i] == ready_df.sacc[i + 1]:
            to_sum.append(ready_df.gene_length_and_intergenic_space[i])
        else:
            staxids.append(ready_df.staxid[i])
            saccs.append(ready_df.sacc[i])
            to_sum.append(ready_df.gene_length_and_intergenic_space[i])
            product = 0
            for number in to_sum:
                product = product + number
            total_space.append(product)
            to_sum = []
    df_of_total_space = pd.DataFrame(list(zip(staxids, saccs, total_space)), columns=['staxid', 'sacc', 'total_space'])
    # total space calculated for each accession
    # new dataframe is created using staxid, sacc and total space
    saccs = []
    column = []
    count = 2
    for i in range(len(df_of_total_space.staxid)):
        if i == 0:
            column.append(1)
            saccs.append(df_of_total_space.sacc[i])
        elif df_of_total_space.staxid[i] == df_of_total_space.staxid[i - 1]:
            if df_of_total_space.sacc[i] in saccs:
                continue
            else:
                column.append(count)
                saccs.append(df_of_total_space.sacc[i])
                count += 1
        else:
            column.append(1)
            saccs.append(df_of_total_space.sacc[i])
            count = 2
    df_of_total_space['sacc_count'] = column
    df_of_total_space.to_csv(r'/home/centos/project/dfs/{}_total_space.csv'.format(name))
    # new metric created - chromosome count for each strain
    # total_space dataframe is exported as .csv

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
    genes = ['pul1', 'pul2', 'pul3', 'pul4']
    grouped = df.groupby(df.qacc)
    pul1 = grouped.get_group("pul1")
    pul2 = grouped.get_group("pul2")
    pul3 = grouped.get_group("pul3")
    pul4 = grouped.get_group("pul4")
    pul1.to_csv(r'/home/centos/project/dfs/pul1.csv')
    pul2.to_csv(r'/home/centos/project/dfs/pul2.csv')
    pul3.to_csv(r'/home/centos/project/dfs/pul3.csv')
    pul4.to_csv(r'/home/centos/project/dfs/pul4.csv')
    # master df split by gene
    cluster_df_names = []
    for gene in genes:
        target_frame = locals()[gene]
        grouped = target_frame.groupby(target_frame.clusters)
        names_of_unique_clusters = target_frame.clusters.unique()
        unwanted = ['qacc', 'pident', 'qcovs', 'sstart', 'send', 'qseq', 'clusters']
        for item in unwanted:
            del target_frame[item]
        # unwanted columns are ditched
        for name in names_of_unique_clusters:
            df_name = gene + '_' + 'cluster' + name
            cluster_df_names.append(df_name)
            locals()[df_name] = grouped.get_group(name)
            # make string into variable name
            # gene dataframes split by cluster
    count = 0
    for name in cluster_df_names:
        target_frame2 = locals()[name]
        target_frame2 = target_frame2.reset_index(drop=True)
        if len(target_frame2.staxid) > 1:
            DNDS()
            target_frame2.to_csv(r'/home/centos/project/dfs/dnds_{}.csv'.format(name))
        else:
            with open(r'/home/centos/project/dfs/error_statement{}.txt'.format(count), 'w') as f:
                f.write('This({}) dataframe has only has one entry'.format(name))
                count += 1
    # DNDS calculated for each dataframe split on gene and cluster
    clusters = []
    for i in range(prediction):
        clusters.append(str(i))
    grouped = df.groupby(df.clusters)
    cluster_df_names2 = []
    for cluster in clusters:
        name = 'cluster{}'.format(cluster)
        cluster_df_names2.append(name)
        locals()[name] = grouped.get_group(cluster)
    # master dataframe split on just clusters
    for name in cluster_df_names2:
        ready_df = locals()[name]
        unwanted = ['pident', 'qcovs', 'qseq', 'sseq','clusters']
        for item in unwanted:
            del ready_df[item]
        intergenic_space_count()
        GenePair()
        total_space()
        culprits = ready_df[ready_df['gene_pair'] == 'NA'].index
        ready_df.drop(culprits, inplace = True)
        if len(ready_df.qacc) > 0:
            ready_df.to_csv(r'/home/centos/project/dfs/{}.csv'.format(name))
    # intergenic sequence length calculated for each dataframe split by cluster
    # each intergenic sequence is assigned a gene pair
    R_script('/home/centos/project/code/Script1.R')
    # R script is run


    # account_for_Ns()
    # df.to_csv(r'/home/centos/project/dfs/df1.csv')
    # R_script('/home/centos/project/code/Script1.R')
    # batch_from_string()
    # database_from_string()
    # make_blast_db('custom2.txt', 'custom2', 'tax_map2.txt')
    # batch_blastn('batch2.txt', 'custom2', 'output2.txt')
    # dataframe(r'/home/centos/blast+/outputs/output2.txt')
    # hist_me_up(r'/home/centos/project/outputs/histogram3.png')
    # if empty_df == False:
    #     R_script('/home/centos/project/code/Script2.R')
    # calling functions
    # data frames are exported as .csv files
    # Second R script is run if data frame contains data
