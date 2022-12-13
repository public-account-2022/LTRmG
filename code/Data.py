import random
import pandas as pd
import operator
import numpy as np


def DownSample(n, df):
    '''sample only part of cells from all sequenced cells'''

    cell_names = list(set(df.loc[:,"cellID"].tolist()))
    if len(cell_names) <= n:
        return range(0,df.shape[0])
    index = range(0,len(cell_names))
    s = random.sample(index,n)
    cell_names_array = np.array(cell_names)
    sample = cell_names_array[s].tolist()
    sample_index =[]
    for i in range(0,df.shape[0]):
        if df.iloc[i,0] in sample:
            sample_index.append(i)
    return sample_index

# take all the cells from PT, Lung and sample 200 cells from each of the other tissues as a representation for determining clones
def ReadFile(inputfile,sample)->(pd.DataFrame, list, pd.DataFrame):
    '''import, filter, clean input file to ultimately generate a distance matrix used for clustering'''

    sample_size = 200
    tissues = sample.loc[:,'tissues'].tolist()
    b = sample.loc[:,'sample'].tolist()
    for i in range(len(b)):
        if b[i] == 't':
            b[i] = True
        else:
            b[i] = False
    precluster_matrix = pd.DataFrame(columns=(['cellID', 'intID','sample']))
    all_data = inputfile.loc[:,['cellID', 'intID','sample']]
    for i in range(len(b)):
        if b[i]:
            tissue = tissues[i]
            to_sample = all_data.loc[all_data['sample']==tissue,:]
            sampled = to_sample.iloc[DownSample(sample_size,to_sample),:]
            precluster_matrix = pd.concat([precluster_matrix,sampled])
        else:
            tissue = tissues[i]
            sampled = all_data.loc[all_data['sample'] == tissue, :]
            precluster_matrix = pd.concat([precluster_matrix, sampled])
    #print(precluster_matrix)
    precluster_matrix = precluster_matrix.drop_duplicates(['cellID', 'intID'], keep='first')
    precluster_matrix.loc[:, 'presence'] = 1
    precluster_matrix = precluster_matrix.loc[:, ['cellID', 'intID', 'presence']]
    print("number of cell-barcode pairs ", precluster_matrix.shape[0])
    # built cell-barcode table, which demonstrate what combination of barcodes each cell have
    cells = list(set(precluster_matrix.loc[:, "cellID"].tolist()))
    print("number of cells ", len(cells))
    static_barcodes = list(set(precluster_matrix.loc[:, "intID"].tolist()))
    print("number of unique static barcodes ", len(static_barcodes))
    cell_barcodes = pd.DataFrame(index=static_barcodes, columns=cells)
    cell_barcodes.fillna(0, inplace=True)
    for i in range(precluster_matrix.shape[0]):
        cell_barcodes.loc[precluster_matrix.iloc[i, 1], precluster_matrix.iloc[i, 0]] = 1
    # put the cells with the exact same combination of barcodes into the same pre-cluster
    pre_cluster = []
    for i in cells:
        barcode = cell_barcodes.loc[:, i].tolist()
        a = True
        k = 0
        while k < len(pre_cluster) and a:
            rep = cell_barcodes.loc[:, pre_cluster[k][0]].tolist()
            if operator.eq(rep, barcode):
                a = False
                pre_cluster[k].append(i)
            k = k + 1
        if a:
            pre_cluster.append([i])
    print("identical static barcode clusters ", len(pre_cluster))
    # calculate the number of overlapped barcode between two cell clusters (time needed)
    print("generating distance matrix")
    cluster_matrix = pd.DataFrame(index=range(len(pre_cluster)), columns=range(len(pre_cluster)))
    for i in range(len(pre_cluster)):
        for x in range(len(pre_cluster)):
            cluster_matrix.iloc[i, x] = OverLap(cell_barcodes.loc[:, pre_cluster[i][0]].tolist(),
                                                cell_barcodes.loc[:, pre_cluster[x][0]].tolist())
    # calculate the distance between each pair of cell clusters and put it in a new maxtrix (time needed)
    cluster_matrix_distance = pd.DataFrame(index=range(len(pre_cluster)), columns=range(len(pre_cluster)))
    for i in range(len(pre_cluster)):
        for x in range(len(pre_cluster)):
            cluster_matrix_distance.iloc[i, x] = (1 - 2 * (
                        cluster_matrix.iloc[i, x] / (cluster_matrix.iloc[i, i] + cluster_matrix.iloc[x, x]))) * 100
    return cluster_matrix_distance, cells, pre_cluster

def OverLap(a, b)-> int:
    '''calculate the number of identical sattic barcode between two cells'''
    sum = 0
    for i in range(len(a)):
        sum = sum + a[i]*b[i]
    return sum

