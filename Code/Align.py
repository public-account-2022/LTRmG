# %%
# load the libraries 
import numpy as np
import numpy.typing as npt
import pandas as pd
import networkx as nx
from threading import Thread
import datetime

# %%
def align(reference: npt.NDArray, barcodes: npt.NDArray, MATCH:int, MISMATCH:int, INDEL:int) -> npt.NDArray:
    '''
    The function called to calculate the alignment score for two evolving barcodes.
    '''
    n = len(reference)
    # one score array for each barcode
    scores = []
    for i in range(n):
        scores.append([])

    try:
        # try to use multiple threads mode
        print("In multiple threads mode.")
        for i in range(n):
            Thread(target = Score, args = (reference[i], barcodes[i], scores[i], MATCH, MISMATCH, INDEL)).run()
    except Exception as e:
        # use single thread mode
        print(e)
        print("Error: unable to start threads. \n Restart in single thread mode.")
        scores = []
        for i in range(n):
            scores.append([])
        for i in range(n):
            Score(reference[i], barcodes[i], scores[i], MATCH, MISMATCH, INDEL)
        
    return np.asarray(scores).reshape(n)
    

# %%
def Score(reference: str, barcode: str, result, MATCH:int, MISMATCH:int, INDEL:int) -> None:
    '''
    The function called to calculate the alignment score for two sequences.
    '''
    refLen = len(reference)
    barLen = len(barcode)

    # the default match, mismatch, insertion or deletion score are as follows
    # MATCH = 1
    # MISMATCH = -2
    # INDEL = -1
    # however, users can specify these parameters in the command line

    # calculate the lowest score (cannot reach actually)
    maxPenal = min(MATCH, MISMATCH, INDEL) 
    maxLen = max(refLen, barLen)

    # record the maximum score and the way to it
    score = np.full((barLen + 1, refLen + 1), maxPenal * maxLen, dtype = int) 
    how = np.full((barLen + 1, refLen + 1), {})

    # initialize the score when if gap
    score[0][0] = 0 
    for i in range(1, barLen + 1):
        score[i][0] = INDEL * i
    for i in range(1, refLen + 1):
        score[0][i] = INDEL * i

    # dynamic programming with Needleman-Wunsch algorithm for global alignment
    for i in range(1, barLen + 1):
        for j in range(1, refLen + 1):
            # if the corresponding position of barcode equals reference
            if (barcode[i-1] == reference[j-1]):
                score[i][j] = max(score[i][j], score[i-1][j-1] + MATCH)
            # if the corresponding position of barcode does not equal reference
            else:
                score[i][j] = max(score[i][j] + MISMATCH, score[i-1][j] + INDEL, score[i][j-1] + INDEL)
    
    # return the maximum alignment score
    result.append(score[barLen, refLen])


#%%