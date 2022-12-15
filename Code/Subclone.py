# %%
# load the packages
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
import random
import Align as an
import Clustering as cu
import Visual as vi
import circlify as circ


# %%
def Seperate_bycloneID(load: pd.DataFrame, data: npt.NDArray) -> dict:
    '''
    The function called to seperate cells by cloneID,
    and attach cellID with intID (static barcodes) and evolving barcodes.
    '''
    cloneID = load.cloneID.unique()  # extract out cloneID
    cloneID_cellID = {key: set() for key in cloneID}
    cellID = load.cellID.unique()  # extract out cellID
    cellID_intID_seq = {key: [] for key in cellID}

    for i in range(len(data)):
        cloneID = data[i][2]  # scan the cloneID
        cellID = data[i][1]  # scan the cellID
        intID = data[i][0]  # scan the intID
        sequences = data[i][3:8]  # scan the sequence1-5
        intID_sequences = {intID: sequences}
        cloneID_cellID[cloneID].add(cellID)
        cellID_intID_seq[cellID].append(intID_sequences)

    return (cloneID_cellID, cellID_intID_seq)


# %%
def Choose_randombarcodes(cloneID_cellID: dict, cellID_intID_seq: dict, random_seed = 200) -> dict:
    '''
    The function called to choose one random barcode from each clone as the reference for barcode alignment.
    '''    
    # record barcodes for each clone
    random_barcodes = {}
    # choose barcode by random 
    for i in cloneID_cellID.keys():
        clone_i_cellID = cloneID_cellID[i].copy()
        random_barcode_cellIDs = list(clone_i_cellID)
        random_barcode_cellIDs.sort()  # sort will ensure the reproducible from set to list 
        # choose the random cellID
        random.seed(random_seed)
        random_barcode_cellID = random_barcode_cellIDs[random.randint(0, len(clone_i_cellID)-1)]
        # choose the random intID 
        random.seed(random_seed)
        random_barcode_intID = cellID_intID_seq[random_barcode_cellID][random.randint(0, len(cellID_intID_seq[random_barcode_cellID])-1)]
        # record the random barcode for each clone 
        random_barcode = list(random_barcode_intID.values())[0]
        random_barcodes[i] = random_barcode
    
    return random_barcodes


# %%
def Barcode_Alignment(cloneID_cellID: dict, cellID_intID_seq: dict, random_barcodes: dict,MATCH: int,MISMATCH:int,INDEL: int) -> dict:
    '''
    The function called to align evolving barcodes of each cells to the reference barcode by each clone.
    '''
    # extract out cellID to prepare a recording list 
    cellIDs = []
    for cloneID in cloneID_cellID.keys():
        cellID = list(cloneID_cellID[cloneID])
        cellIDs.extend(cellID)
    cellID_intID_score={key: [] for key in cellIDs}

    # scoring each barcode with the reference barcode in each cell for each clone 
    for cloneID in cloneID_cellID.keys():
        cellID_i=list(cloneID_cellID[cloneID])
        for cellID in cellID_i:
            intID_seq=cellID_intID_seq[cellID]
            for i in range(len(intID_seq)):
                intID=list(intID_seq[i].keys())[0]
                seq=list(intID_seq[i].values())[0]
                score=an.align(reference = random_barcodes[cloneID], barcodes = seq, MATCH=MATCH,MISMATCH=MISMATCH,INDEL=INDEL)
                intID_score={intID: score}
                cellID_intID_score[cellID].append(intID_score)

    return cellID_intID_score


# %%
def Cell_Clustering(cloneID_cellID: dict, cellID_intID_score: dict) -> dict:
    # record subclones for each clone 
    subclones = {}
    for cloneID in cloneID_cellID.keys():
        cellID_i = list(cloneID_cellID[cloneID])  # extract out cellIDs corresponding to each cell
        sub_cellID_intID_score = {}
        for cellID in cellID_i:
            sub_cellID_intID_score[cellID] = cellID_intID_score[cellID]
        subclone = cu.AGNES(sub_cellID_intID_score)  # clustering for each clone 
        subclones[cloneID] = subclone
    
    return subclones


# %%
 