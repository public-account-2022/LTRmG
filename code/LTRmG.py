# load the basic packages
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import Subclone as sc
import Visual as vi
import sys
import time
import pandas as pd
from CloneDefination import HierarchicalClustering
from Data import ReadFile
from optparse import OptionParser
from CloneFilter import FilteredClone
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", default='example.txt',
                  help="write report to FILE", metavar="FILE")
parser.add_option("-o", "--output", dest="outputpath", default = './results/subclones.txt',
                  help="specify the output path", metavar="OUTPUT")
parser.add_option("-s", "--sample", dest="sample", default = 'samples.txt',
                  help="a table which specify which tissue needs to be sample", metavar="SAMPLE")
parser.add_option("-n", "--numberofclusters", dest="numberofclusters", default = 68,
                  help="the initial amount of clusters, which depends on the experiment and dataset, should be bigger than hypothetical amount of clusters for over-clustering", metavar="NUMBEROFCLUSTERS")
parser.add_option("-a", "--match", dest="MATCH", default = 1, help = 'score for each match',metavar="MATCH")
parser.add_option("-b", "--mismatch", dest="MISMATCH", default = -2,help = 'score for each mismatch', metavar="MISMATCH")
parser.add_option("-c", "--indel", dest="INDEL", default = -1,help = 'score for each indel', metavar="INDEL")
parser.add_option("-p", "--position", dest="primaryposition", default = 'PT',help = 'state which sample is the original site of cells', metavar="POSITION")
#you need to manually set the variable 'parameters', if you want to run this program in python, following a structure like['-f','inputfile.txt',……]
parameters = sys.argv[1:]

#load file from command line input
(options, args) = parser.parse_args(parameters)
filename = options.filename
outputpath = options.outputpath
numberofclusters = int(options.numberofclusters)
samplefile = options.sample
MATCH = int(options.MATCH)
MISMATCH = int(options.MISMATCH)
INDEL = int(options.INDEL)
primary_position = options.primaryposition
inputfile =pd.read_table(filename,sep = '\t',header = 0)
begin_time = time.time()
if 'cloneID' in inputfile.columns:
    singlet = inputfile
else:
    sample = pd.read_table(samplefile, sep='\t', header=0)
    begin_time1 = time.time()

    # data process for cluster identification of clones
    cluster_matrix_distance, cells, pre_cluster = ReadFile(inputfile, sample)

    # identify cell clone based on Hierarchical Cluster
    cell_clone = HierarchicalClustering(cluster_matrix_distance, numberofclusters, cells, pre_cluster)

    # filter the clone and select singlet
    singlet = FilteredClone(inputfile, cell_clone, primary_position)
    end_time1 = time.time()
    print("clustering time consumption: ", end_time1 - begin_time1, 's')

    singlet["cloneID"] = list(map(int, list(singlet["cloneID"])))

# %%
# For this Lineage Tree Recontruction Python package quick running
# %%
# load the barcodes data
# extract out related columns
load = singlet[['intID', "cellID", "cloneID", "sequence1", "sequence2",
             "sequence3", "sequence4", "sequence5"]]  # pandas.DataFrame

data = load.to_numpy()  # numpy.ndarray
# %%
### Cell Clone Clustering

# Seperate cells to clones by static barcodes (intID)


# %%
### Cell Subclone Clustering

# Seperate cells by cloneID
# Attach cellID with intID (static barcodes) and evolving barcodes
cloneID_cellID, cellID_intID_seq = sc.Seperate_bycloneID(load, data)

# Choose one random barcode from each clone as the reference for alignment
random_barcodes = sc.Choose_randombarcodes(cloneID_cellID, cellID_intID_seq, random_seed=200)

# Align evolving barcodes of each cells to the reference barcodes
print("Now barcodes are aligning, and it may take a few minutes to complete...")
cellID_intID_score = sc.Barcode_Alignment(cloneID_cellID, cellID_intID_seq, random_barcodes,MATCH,MISMATCH,INDEL)

# Calculate the average distance and do AGNES (AGglomerative NESting) clustering for those cells
# The distance between cells is defined as the number of barcodes + intID comparison + alignment score
subclones = sc.Cell_Clustering(cloneID_cellID, cellID_intID_score)


# %%
### Visualization of Clustering Results

# Choice 1 Circle Plot

# You need to prepare the data for circle plot visualization
circle = vi.prepare_circle(load, subclones)

# output the clone and subclone results
os.makedirs("./results")
f = open(outputpath, "w")
f.writelines(str(circle))
f.close()

# View clone information only
vi.circle_allclones(circle)

# View clone and sublone information
vi.circle_allclones_subclones(circle)
end_time = time.time()
print('total time consumption ', end_time-begin_time,'s')

# View subclone information for only one clone
# vi.circle_single_subclone(subclones, cloneID=0)

# %%
# Chocie 2 Edge Lists

## You can choose the align argument: 'vertical' or 'horizontal'


# View edge lists for each clone
#vi.Edgelists_Clone(subclones, align="vertical")

# View edge list for only one clone
#vi.Edgelist_OneClone(cloneID=0, subclones=subclones, align="horizontal")

# %%
# Choice 3 Circle Visualization of Edge Lists
#vi.circular_tree(subclones, cloneID=0)

# %%
# Choice 4 Lineage Tree Visualization
#vi.node_tree(subclones, cloneID=0)






