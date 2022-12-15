# %%
# load the packages
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
import queue
from collections import Counter


# %%
# Calculate the distance for clustering 

def Barcodes_Distance(cell_a_cellID: str, cell_b_cellID: str, sub_cellID_intID_score: dict) -> float:
    '''
    The function called to calculate and return a distance between two cells from different clusters. 
    '''
    # define the scroe recording matrix 
    record_matrix = [0 for i in range(1)]

    # dealing with data structures 
    cell_a = sub_cellID_intID_score[cell_a_cellID] 
    cell_b = sub_cellID_intID_score[cell_b_cellID]

    cell_a_copy = np.asarray([[list(i.keys())[0], list(i.values())[0]] for i in cell_a], dtype=object)
    cell_a_copy = pd.DataFrame(cell_a_copy).sort_values(by=0).to_numpy()  # sort intIDs for comparison
    cell_b_copy = np.asarray([[list(i.keys())[0], list(i.values())[0]] for i in cell_b], dtype=object)
    cell_b_copy = pd.DataFrame(cell_b_copy).sort_values(by=0).to_numpy()  # sort intIDs for comparison 
    
    # extract out cell related intID data 
    cell_a_intID = [list(i.keys())[0] for i in cell_a] 
    cell_b_intID = [list(i.keys())[0] for i in cell_b]

    cell_a_intID_number = len(cell_a_intID)
    cell_b_intID_number = len(cell_b_intID)

    # the difference in the number of barcodes 
    record_matrix[0] = abs(cell_a_intID_number - cell_b_intID_number)  

    # different situations tha decide cell-cell distance 
    if cell_a_intID_number == cell_b_intID_number:  # same intID number 
        a = Counter(cell_a_intID)
        b = Counter(cell_b_intID)

        if a == b:  # same intID 
            record_matrix[0] = 0  # no different intID 

            # calculate differences between alignment scores 
            differences = []  
            for i in range(len(cell_a_copy)):
                difference = abs(cell_a_copy[i][1] - cell_b_copy[i][1])
                differences.append(difference)
            record_matrix.extend(np.asarray(differences).mean(axis=0))

            distance = 5*record_matrix[0] + sum([i*1 for i in record_matrix[1:]])  # 0.5 + 0.1 * 5 
            return distance

        elif a != b:  # different intID
            inter = set(cell_a_intID).intersection(set(cell_b_intID))  # take the intersection of two intID sets 
            
            if inter == set():  # no intID in common 
                record_matrix[0] = cell_a_intID_number  # all are different 
                record_matrix.extend(np.asarray([0 for i in range(5)]))  # align score = 0
                distance = 5*record_matrix[0] + sum([i*1 for i in record_matrix[1:]])  # 0.5 + 0.1 * 5 
                return distance
            else:
                # extract out barcodes with the same intIDs 
                sub_a = np.asarray([i for i in cell_a_copy if i[0] in inter])
                sub_b = np.asarray([i for i in cell_b_copy if i[0] in inter])

                # more different intID
                record_matrix[0] = max(len(cell_a_copy)-len(sub_a), len(cell_b_copy)-len(sub_b))  

                # calculate differences between alignment scores
                differences = []
                i = 0
                count = 0
                # add pointers to make sure that only the barcodes with the same intID will be compared 
                while count <= record_matrix[0] and i < min(len(sub_a), len(sub_b)):
                    if sub_a[i][0] == sub_b[i][0]:
                        count += 1
                        difference = abs(sub_a[i][1] - sub_b[i][1])
                        differences.append(difference)
                        i += 1
                    else:
                        i += 1                    
                record_matrix.extend(np.asarray(differences).mean(axis=0))

                distance = 5*record_matrix[0] + sum([i*1 for i in record_matrix[1:]])  # 0.5 + 0.1 * 5 
                return distance

    elif cell_a_intID_number != cell_b_intID_number:  # different intID number 
        inter = set(cell_a_intID).intersection(set(cell_b_intID))
        if inter == set():  # no intID in common
            record_matrix[0] = max(cell_a_intID_number, cell_a_intID_number)  # max the difference
            record_matrix.extend(np.asarray([0 for i in range(5)]))
            distance = 5*record_matrix[0] + sum([i*1 for i in record_matrix[1:]])  # 0.5 + 0.1 * 5 
            return distance
        else:
            # extract out barcodes with the same intIDs
            sub_a = np.asarray([i for i in cell_a_copy if i[0] in inter])
            sub_b = np.asarray([i for i in cell_b_copy if i[0] in inter])

            # more different intID
            record_matrix[0] = max(len(cell_a_copy)-len(sub_a), len(cell_b_copy)-len(sub_b))

            # calculate differences between alignment scores
            differences = []
            i = 0
            count = 0
            # add pointers to make sure that only the barcodes with the same intID will be compared
            while count <= record_matrix[0] and i < min(len(sub_a), len(sub_b)):
                if sub_a[i][0] == sub_b[i][0]:
                    count += 1
                    difference = abs(sub_a[i][1] - sub_b[i][1])
                    differences.append(difference)
                    i += 1
                else:
                    i += 1
                
            record_matrix.extend(np.asarray(differences).mean(axis=0))

            distance = 5*record_matrix[0] + sum([i*1 for i in record_matrix[1:]])  # 0.5 + 0.1 * 5 
            return distance


def Cell_Distance(clu_cellID_A: list, clu_cellID_B: list, sub_cellID_intID_score: dict) -> npt.NDArray:
    '''
    The function called to compute the average distance between two clusters. 
    It will record the distance of two cells from different clusters. 
    '''
    # extract cellID from two clusters 
    cluster_A = {key: sub_cellID_intID_score[key] for key in clu_cellID_A}
    cluster_B = {key: sub_cellID_intID_score[key] for key in clu_cellID_B}

    all_distance = []
    for cell_a in cluster_A:
        for cell_b in cluster_B:
            cell_distance = Barcodes_Distance(cell_a, cell_b, sub_cellID_intID_score)
            all_distance.append(cell_distance)
    return np.asarray(all_distance).mean(axis=0)  # average distance used between two clusters 


def Average_Distance(record: dict, sub_cellID_intID_score: dict) -> npt.NDArray:
    '''
    The function called to compute the distance matrix between all clusters.
    '''
    # n = the number of clusters at current clustering stage 
    n = len(record)

    if n == 1:  # only one cluster 
        return np.asarray([0])
    else:
        # initialization of the distance matrix 
        dis_matrix = [[0 for i in range(n)] for j in range(n)]

        for key in record.keys():
            clu_cellID_A = record[key]
            for anokey in record.keys():
                if anokey == key:
                    distance = float("inf")  # distance of itself will be infinity for calculation convenience
                else:
                    clu_cellID_B = record[anokey]
                    distance = Cell_Distance(clu_cellID_A, clu_cellID_B, sub_cellID_intID_score)
                dis_matrix[key][anokey] = distance

        return np.asarray(dis_matrix)

# %%
# Implementation of AGglomerative NESting algorithm (AGNES)
def AGNES(sub_cellID_intID_score: dict) -> dict:
    '''
    The function called to cluster these cells into different levels of clusters.
    Based on AGglomerative NESting algorithm (AGNES). 
    This function will return a dictionary that contain clustering results for every iteration. 
    '''
    # initialization (each cell is considered as one cluster)
    cellID = list(sub_cellID_intID_score.keys())
    keys = [i for i in range(len(cellID))]
    # record clustering results for each iteration 
    record = {key: [cellID] for key, cellID in zip(keys, cellID)}  # clu_name: cellID 
    
    num_clu = len(record.keys())  # initialize the original number of clusters 

    if num_clu == 1:  # if only one cell within this clone 
        storage = {0: ({0: cellID}, [0])} # {times: ({clusID: cellID}, [includeLastID])}
        return storage
    
    dis_matrix = Average_Distance(record, sub_cellID_intID_score)
    return dendrogram(dis_matrix, record)

def belong_root(belong, x):
    if (belong[x]==x):
        return x
    belong[x] = belong_root(belong, belong[x])
    return belong[x]


def dendrogram_merge_list(dis_matrix):
    n = len(dis_matrix[0])
    que = queue.PriorityQueue()

    # Init
    dis = []
    for i in range(n*2):
        dis.append([0]*n*2)
    for i in range(n):
        for j in range(n):
            que.put([dis_matrix[i][j],(i,j)])
            dis[i][j]=dis_matrix[i][j]

    new_id = n-1

    flg = [1]*n + [0]*n # whether the point exists (for the PriorityQueue)
    sz = [1]*n + [0]*n # size of each cluster
    ret = [] # return the id of the merged clusters
    # id of newly merged cluster is n, n+1, n+2, etc. 
    # i.e., the id of the cluster created (merged) in time i is n+i-1

    for time in range(1,n):
        tmp = que.get()
        while(not(flg[tmp[1][0]] and flg[tmp[1][1]])): # check whether edge is still valid
            tmp = que.get()

        # merge cluster x, y
        x,y = tmp[1]
        new_id += 1
        sz[new_id] = sz[x] + sz[y]
        flg[x] = flg[y] = 0
        flg[new_id] = 1
        ret.append((x,y))

        # update distance between the new point and others
        for j in range(new_id):
            if(flg[j]):
                dis[new_id][j] = dis[x][j] + dis[y][j]
                dis[j][new_id] = dis[new_id][j]
                que.put([dis[new_id][j]/sz[new_id]/sz[j], (new_id, j)])

    return ret

def dendrogram_build_result(merge_list, record):
    n = len(record)
    now = {}
    for i in range(n):
        now[i] = record[i] 
    result = {0: (now, [[i] for i in range(n)])}
    new_id = n-1
    belong = [i for i in range(2*n)]
    for time in range(len(merge_list)):
        x,y = merge_list[time]
        x = belong[x]
        y = belong[y]
        new_id += 1
        new = {}
        merge = [[x, y]]
        new[0] = now[x] + now[y]
        num = 1
        mapp = {}
        for key in now.keys():
            if (not(key == x or key == y)):
                mapp[key] = num
                new[num] = now[key]
                merge.append([key])
                num += 1

        for i in range(len(belong)):
            if(belong[i] in mapp):
                belong[i] = mapp[belong[i]]

        belong[new_id] = 0

        now = new
        result[time+1] = (now, merge)
    return result


def dendrogram(dis_matrix, record):
    merge_list = dendrogram_merge_list(dis_matrix)
    return dendrogram_build_result(merge_list, record)


# %%
# Implementation of AGglomerative NESting algorithm (AGNES)

def AGNES_old(sub_cellID_intID_score: dict) -> dict:
    '''
    The function called to cluster these cells into different levels of clusters.
    Based on AGglomerative NESting algorithm (AGNES). 
    This function will return a dictionary that contain clustering results for every iteration. 
    '''
    # initialization (each cell is considered as one cluster)
    cellID = list(sub_cellID_intID_score.keys())
    keys = [i for i in range(len(cellID))]
    # record clustering results for each iteration 
    record = {key: [cellID] for key, cellID in zip(keys, cellID)}  # clu_name: cellID 
    
    cur_clu = len(record.keys())  # initialize the original number of clusters 

    if cur_clu == 1:  # if only one cell within this clone 
        storage = {0: ({0: cellID}, [0])}
        return storage

    else:
        times = 0  # iteration times

        # store the results of each time of clustering 
        storage = {0: (record, [[i] for i in range(len(cellID))])}  
        
        while cur_clu > 1:  # stop when become the largest cluster 

            # compute the average distance
            dis_matrix = Average_Distance(record, sub_cellID_intID_score)

            # update the assignment 
            min_dis = dis_matrix.min(axis=1)

            # link this clsuter to its most closed related cluster(s)
            operation = {}  # clu_name: minimum clu_name 
            for cluster in range(len(min_dis)):
                operation[cluster] = [i for i in range(len(dis_matrix[cluster])) if dis_matrix[cluster][i] == min_dis[cluster]]

            # form new clusters based on previous clusters 
            clusters = []  
            counter = {i: 0 for i in range(len(record))}  # each cluster only appear once 
            for clu_name in operation.keys():
                if counter[clu_name] > 0:  # no repeat
                    continue
                else:
                    next_clu = [clu_name]
                    counter[clu_name] += 1
                    corr_clu_name = operation[clu_name]
                    for clu in corr_clu_name:
                        if counter[clu] > 0:  # no repeat 
                            continue
                        else:
                            next_clu.append(clu)
                            counter[clu] += 1
                    clusters.append(next_clu)
            
            # assign cellID to each new cluster
            values = []
            for new_clu_name in range(len(clusters)): 
                cluster = clusters[new_clu_name]
                cluster_cellID = []
                for key in cluster:
                    cluster_cellID_list = storage[times][0][key]
                    cluster_cellID.extend(cluster_cellID_list)
                # cluster_cellID = [storage[times][key] for key in cluster] but nested lists 
                values.append(cluster_cellID)
            keys = [i for i in range(len(clusters))]
            record = {key: value for key, value in zip(keys, values)}
            
            
            # update iteration times
            times += 1

            # record clustering results for each iteration 
            storage[times] = (record, clusters)

            # update the number of clusters for the judgement 
            cur_clu = len(record.keys())
    
        return storage
