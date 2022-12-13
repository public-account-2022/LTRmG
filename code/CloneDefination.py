import random
import pandas as pd
import operator
import numpy as np
def Clusterdistance(df,a,b)->float:
    ''' calculate the distance between two clusters with self-developed method'''

    maxvaluesA=[]
    maxvaluesB=[]
    for n in a:
        d = 0
        for m in b:
            if df.iloc[n,m] >d :
                d = df.iloc[n,m]
        maxvaluesA.append(d)
    meanmaxA=sum(maxvaluesA)/len(maxvaluesA)
    for m in b:
        d = 0
        for n in a:
            if df.iloc[n,m] >d :
                d = df.iloc[n,m]
        maxvaluesB.append(d)
    meanmaxB = sum(maxvaluesB) / len(maxvaluesB)
    result=(meanmaxA+meanmaxB)/2
    return result

def HierarchicalClustering(df, k, cells, pre_cluster)->pd.DataFrame:
    '''cluster the identical static barcode clusters(cells with identical barcode combination were clustered together in advance) with Hierarchical Clustering method'''

    print("clustering")
    clusters = []
    origin_distance = df.copy()
    for i in range (df.shape[0]):
        clusters.append([i])
    while len(clusters) > k: #each round two clusters with the smallest inter-cluster distance are merged
        distance = 1000
        to_merge = []
        for x in range(len(clusters)):
            for y in range(len((clusters))):
                if x !=y :
                    cd = df.loc[x,y]
                    if cd < distance:
                        distance = cd
                        to_merge = [x,y]
        clusters[to_merge[0]].extend(clusters[to_merge[1]]) #rebuild the clusters list and distance matrix after cluster collapse
        del clusters[to_merge[1]]
        df = df.drop(columns=[to_merge[1]],axis=1)
        df = df.drop(index=[to_merge[1]], axis=0)
        df.index = range(0,len(clusters))
        df.columns = range(0,len(clusters))
        for m in range(0,len(clusters)):
            df.iloc[m,to_merge[0]] = Clusterdistance(origin_distance,clusters[m],clusters[to_merge[0]])
            df.iloc[to_merge[0],m] = Clusterdistance(origin_distance,clusters[to_merge[0]],clusters[m])
    label = [0]*origin_distance.shape[0]
    for m in range(0,len(clusters)):
        index = clusters[m]
        for n in index:
            label[n] = m
    print("cluster finished")
    cell_clone = pd.DataFrame(columns=(["cellID", "cloneID"]))
    cell_clone.loc[:, "cellID"] = cells
    for i in range(len(pre_cluster)): # assign each cell to one clone
        cell_same_barcode = pre_cluster[i]
        for x in cell_same_barcode:
            cell_clone.loc[cell_clone['cellID'] == x, "cloneID"] = label[i]
    print(cell_clone)
    return cell_clone