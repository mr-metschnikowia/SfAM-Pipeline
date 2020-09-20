import os

folders = [r'/home/centos/project/outputs/',r'/home/centos/project/needle_output/',r'/home/centos/project/genomes/',r'/home/centos/project/dfs/'
           ,r'/home/centos/blast+/outputs/',r'/home/centos/blast+/dbs/',r'/home/centos/blast+/outputs/blastn_output/',r'/home/centos/project/']
for folder in folders:
    for file in os.listdir(folder):
            path = folder + file
            try:
                os.remove(path)
            except IsADirectoryError:
                continue
