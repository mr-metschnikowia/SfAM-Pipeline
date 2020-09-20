# pairwise global alignment occurs between pul gene DNA from each strain and apc1.2 pul genes
# > data from all alignments is processed and returned as master table in following format: qacc,staxid,sacc,pident,qcovs,sstart,send,qseq,sseq

import os
import pandas as pd
# importing prereqs

def needle(species):
    os.system('needleall -asequence /home/centos/project/genes/batch.txt -bsequence /home/centos/project/genomes/{}_genes.txt -outfile /home/centos/project/needle_output/{}.csv -gapopen 10 -gapextend 0.5'.format(species,species))
    # global batch pairwise alignment

def csv_convert():
    staxid = 0
    list_of_staxids = []
    list_of_species = []
    for file in os.listdir(r'/home/centos/project/needle_output'):
        list_of_staxids.append(str(staxid))
        species = file[:file.find('.csv')]
        list_of_species.append(species)
        with open(r'/home/centos/project/needle_output/{}'.format(file), 'r') as f:
            data = f.readlines()
            data = data[:len(data) - 3]
        new_lines = []
        for line in data:
            split_line = line.split(' ')
            split_line[3] = split_line[3][1:len(split_line[3]) - 2]
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
        staxid += 1
    staxid_to_species = pd.DataFrame(list(zip(list_of_staxids, list_of_species)), columns=['staxid', 'species'])
    staxid_to_species.to_csv('/home/centos/project/taxmap.csv', index=False)
    # create tax_map

def process_needle_output():
    for file in os.listdir(r'/home/centos/project/needle_output'):
        species = file[:file.find('.csv')]
        # species and staxid are defined
        df = pd.read_csv(r'/home/centos/project/needle_output/{}'.format(file))
        df = df.loc[df.groupby('qacc')['qcovs'].idxmax()]
        # needle output is read for each strain
        # output is filtered by top hits for each gene
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
            sstart.append(int(gene[0]))
            send.append(int(gene[1]))
        df['sstart'] = sstart
        df['send'] = send
        qseq = []
        for i in df.send:
            qseq.append('dummy')
        df['qseq'] = qseq
        # add dummy qseq column
        df['sseq'] = sseq
        df.to_csv(r'/home/centos/project/needle_output/{}'.format(file),index=False)
        # needle output modified

def prepare_protein_output():
    for file in os.listdir('/home/centos/blast+/outputs/'):
        try:
            with open(r'/home/centos/blast+/outputs/{}'.format(file), 'r') as f:
                data = f.readlines()
                new_data = ['qacc\n']
                for line in data:
                    new_line = line[:4]
                    new_data.append(new_line + '\n')
            with open(r'/home/centos/blast+/outputs/{}'.format(file), 'w') as f:
                f.writelines(new_data)
            # title line is inserted
        except IsADirectoryError:
            continue

def remove_trash():
    for file in os.listdir('/home/centos/project/needle_output/'):
        species = file[:file.find('.csv')]
        not_these_ones = ['ucd127','ap47','277','bath1','mp5']
        if species not in not_these_ones:
            path = '/home/centos/project/needle_output/' + file
            df = pd.read_csv(path)
            df_proteins = pd.read_csv('/home/centos/blast+/outputs/{}_proteins.txt'.format(species),usecols=[0])
            quacks = []
            for qacc in df_proteins.qacc:
                quacks.append(qacc)
            for qacc in df.qacc:
                if qacc not in quacks:
                    culprits = df[df['qacc'] == qacc].index
                    df.drop(culprits, inplace=True)
            df.to_csv('/home/centos/project/needle_output/{}.csv'.format(species),index=False)
    # remove intruders

def concatenate():
    count = 0
    for file in os.listdir('/home/centos/project/needle_output/'):
        path = '/home/centos/project/needle_output/' + file
        df = pd.read_csv(path)
        if count == 0:
            switch = True
        else:
            switch = False
        df.to_csv('/home/centos/project/output.csv', mode='a', header=switch, index=False)
        count += 1

def correct_accessions():
    df = pd.read_csv('/home/centos/project/output.csv')
    for i in range(len(df.sacc)):
        df.sacc[i] = df.sacc[i].split('.')[0]
    # removing unique identifier from accession code
    df.to_csv('/home/centos/project/output.csv',index=False)

if __name__ == '__main__':
    needle('ucd127')
    needle('mp5')
    needle('bath1')
    needle('ap47')
    needle('277')
    needle('apc1.2')
    needle('c.auris')
    needle('k.lactis')
    needle('k.dobzhanskii')
    needle('z.mrakii')
    needle('k.pastoris')
    needle('l.thermotolerans')
    needle('k.marxianus')
    needle('s.eubayanus')
    needle('s.kudriavzevii')
    needle('s.paradoxus')
    needle('s.cerevisiae')
    needle('h.burtonii')
    # call for each species/strain
    csv_convert()
    process_needle_output()
    prepare_protein_output()
    remove_trash()
    concatenate()
    correct_accessions()





