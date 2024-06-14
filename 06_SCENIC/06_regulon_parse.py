##############################BEFORE RUNNING THIS SCRIPT##############################

# sys.argv[1] = location of concatenated regulon files
# sys.argv[2] = the result file including regulon and target-genes

# Recommend to run it on the same environment of pyscenic : $ conda activate pyscenic

# Need to install those modules before running it : $ pip install [name of module]

######################################################################################

import os
import sys
import csv

def read_regulon_file(path):
    regulons = {}
    with open(path, 'r') as f:
        next(f)
        next(f)
        next(f) # pass the first three line since they are the index of each columns
        reader = csv.reader(f, delimiter=",")
        for row in reader:
            regulon = row[0] # Take first column called "TF"
            targets = eval(row[8]) # Take the column called "Enrichment, TargetGenes"
            if regulon not in regulons:
                regulons[regulon] = []
            targets = dict(targets).keys()
            regulons[regulon] = regulons[regulon] | targets
        return(regulons)

def write_new_regulon_file(regulons,output_path):
    f = open(output_path,'w')
    for (key,value) in regulons.items():
        f.write(f'{key}\t{";".join(value)}\n')
    f.close()

if __name__ == "__main__":
    regulon_path = sys.argv[1]
    regulons = read_regulon_file(regulon_path)
    write_new_regulon_file(regulons, sys.argv[2])
