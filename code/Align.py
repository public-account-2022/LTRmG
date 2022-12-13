# %%
# load the libraries
import numpy as np
import numpy.typing as npt
import pandas as pd
import networkx as nx


# %%
def align(reference: npt.NDArray, barcodes: npt.NDArray,MATCH: int,MISMATCH:int,INDEL: int) -> npt.NDArray:
    '''
    The function called to calculate the alignment score for two evolving barcodes.
    '''
    n = len(reference)
    scores = []  # one score array for each barcode
    for i in range(n):
        score = Score(reference=reference[i], barcode=barcodes[i], MATCH=MATCH,MISMATCH=MISMATCH,INDEL=INDEL)
        scores.append(score)
    return np.asarray(scores)


# %%
def Score(reference: str, barcode: str, MATCH: int,MISMATCH:int,INDEL: int) -> int:
    '''
    The function called to calculate the alignment score for two sequences.
    '''
    refLen = len(reference)
    barLen = len(barcode)

    # set the defalut match, mismatch, insertion or deletion score
    MIS = MISMATCH

    # calculate the lowest score (cannot reach actually)
    maxPenal = min(MATCH, MIS, INDEL) 
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
                score[i][j] = max(score[i][j], score[i-1][j] + INDEL, score[i][j-1] + INDEL)
    
    # return the maximum alignment score
    return score[barLen, refLen]


#%%