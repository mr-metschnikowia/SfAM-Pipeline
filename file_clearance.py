import os

folders = [r'/home/centos/project/outputs/',r'/home/centos/project/genomes/',r'/home/centos/project/genes/',r'/home/centos/project/dfs/'
           ,r'/home/centos/blast+/outputs/',r'/home/centos/blast+/tax_maps/',r'/home/centos/blast+/pre_taxmaps/',r'/home/centos/blast+/dbs/']
for folder in folders:
    for file in os.listdir(folder):
        if folder == r'/home/centos/project/genomes/' or folder == r'/home/centos/project/genes/':
            if file == 'custom.txt' or file == 'batch.txt':
                path = folder + file
                os.remove(path)
        elif file == 'reference.txt':
            continue
        else:
            path = folder + file
            os.remove(path)
