# pairwise global alignment occurs between pul gene DNA from each strain and apc1.2 pul genes 
# > data from all alignments is processed and returned as master table in following format: qacc,staxid,sacc,pident,qcovs,sstart,send,qseq,sseq

import os
import pandas as pd
# importing prereqs

def needle(species):
    os.system('needleall -asequence /home/centos/project/genes/batch.txt -bsequence /home/centos/project/genomes/{}_genes.txt -outfile /home/centos/project/needle_output/{}.csv -gapopen 10 -gapextend 0.5'.format(species,species))
    # global batch pairwise alignment
    
def process_needle_output():
    staxid = 0
    list_of_staxids = []
    list_of_species = []
    for file in os.listdir(r'/home/centos/project/needle_output'):
        list_of_staxids.append(staxid)
        species = file[:file.find('.csv')]
        list_of_species.append(species)
        with open(r'/home/centos/project/needle_output/{}'.format(file), 'r') as f:
            data = f.readlines()
            data = data[:len(data)-3]
        new_lines = []
        for line in data:
            split_line = line.split(' ')
            split_line[3] = split_line[3][1:len(split_line[3])-2]
            # remove brackets around score
            split_line.insert(1, str(staxid))
            # insert staxid,
            new_line = ','.join(split_line)
            new_line = new_line + '\n'
            new_lines.append(new_line)
            # replace ' ' with ',' i.e. convert to csv
        with open(r'/home/centos/project/needle_output/{}'.format(file), 'w') as f:
            f.writelines(new_lines)
        with open(r'/home/centos/project/needle_output/{}'.format(file), 'r') as f:
            lines = f.readlines()
        lines.insert(0, 'qacc,staxid,sacc,pident,qcovs\n')
        with open(r'/home/centos/project/needle_output/{}'.format(file), 'w') as f:
            f.writelines(lines)
        # title line is inserted
        df = pd.read_csv(r'/home/centos/project/needle_output/{}'.format(file))
        df = df.loc[df.groupby('qacc')['qcovs'].idxmax()]
        # filter by top hits for each gene
        df['qcovs'] = df['qcovs']/df['pident']
        # score is normalised (dividing my length of alignment)
        with open(r'/home/centos/project/genomes/{}_genes.txt'.format(species),'r') as f:
            data = f.read()
        data = data.split('>')
        # each gene is item in list
        genes = []
        for sacc in df.sacc:
            for gene in data:
                result = gene.find(sacc)
                if result > -1:
                    genes.append(gene)
        # genetic information extracted for corresponding sacc
        sstart= []
        send = []
        sseq = []
        for gene in genes:
            cds = gene[gene.find('[gbkey=CDS]') + len('[gbkey=CDS]'):]
            cds = "".join([i for i in cds if i.isalpha()])
            sseq.append(cds)
            # extract cds
            gene = gene[gene.find('location='):]
            gene = gene[:gene.find(' ')]
            for i in range(len(gene)):
                if gene[i].isnumeric() == True:
                    gene = gene[i:]
                    break
            endings = [',', ')', ']']
            for end in endings:
                if gene.find(end) > -1:
                    gene = gene[:gene.find(end)]
                    gene = gene.split('..')
                    break
            # find sstart and send for each hit
            sstart.append(gene[0])
            send.append(gene[1])
        df['sstart'] = sstart
        df['send'] = send
        qseq = []
        for i in df.send:
            qseq.append('n/a')
        df['qseq'] = qseq
        # add dummy qseq column
        df['sseq'] = sseq
        if staxid == 0:
            switch = True
        else:
            switch = False
        df.to_csv('/home/centos/project/output.csv',mode='a',header=switch,index=False)
        # dataframe added to csv
        staxid += 1
        # change staxid for next needle output file
    staxid_to_species = pd.DataFrame(list(zip(list_of_staxids, list_of_species)), columns=['staxid', 'species'])
    staxid_to_species.to_csv('/home/centos/project/taxmap.csv',index=False)
    # create tax_map

if __name__ == '__main__':
    needle('apc1.2')
    needle('c.auris')
    needle('k.lactis')
    # call for each species/strain
    process_needle_output()





