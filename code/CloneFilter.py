import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re as re
import numpy as np



def Get_num_of_cells_in_each_cluster(df:pd.DataFrame) ->pd.DataFrame:
    '''
    Get a dataframe containing the number of cells in each clusters
    '''
    Clusters = df.drop_duplicates(['cellID'], keep='first').loc[:, ["cloneID"]].values.tolist()
    Clusters_count={}
    for i in range(len(Clusters)):
        if Clusters[i][0] not in Clusters_count:
            Clusters_count[Clusters[i][0]]=1
        else:
            Clusters_count[Clusters[i][0]]=Clusters_count.get(Clusters[i][0])+1
    Clusters_df=pd.DataFrame(pd.Series(Clusters_count),columns=['num_cells_in_cluster']).reset_index().rename(columns={'index':'cloneID'})
    return Clusters_df



def Get_num_cells_per_intID_in_cluster(df:pd.DataFrame) ->pd.DataFrame:
    '''
    Get a dataframe containing the num_cells_per_intID_in_cluster
    '''
    DF_CellIntCluster = (df.loc[:, ["cellID", "intID", "cloneID"]]).iloc[:, [1, 2]]
    tmplist = []
    for i in range(DF_CellIntCluster.shape[0]):
        str1 = str(DF_CellIntCluster.iloc[i, 0]) + "_" + str(DF_CellIntCluster.iloc[i, 1])
        tmplist.append(str1)
    DF_CellIntCluster["num_cells_per_intID_in_cluster"] = tmplist
    tmpdf = pd.DataFrame(DF_CellIntCluster.loc[:, "num_cells_per_intID_in_cluster"].value_counts(),
                         columns=["num_cells_per_intID_in_cluster"]).reset_index().rename(
        columns={"index": "ComplexID"})
    intIDs = []
    cloneIDs = []
    for i in range(tmpdf.shape[0]):
        tmprow = list((tmpdf.iloc[i, 0]).split("_"))
        intIDs.append(tmprow[0])
        cloneIDs.append(int(float(tmprow[1])))
    tmpdf["intID"] = intIDs
    tmpdf["cloneID"] = (cloneIDs)
    return tmpdf


def Get_fraction_cells_with_intID_in_cluster(df:pd.DataFrame) ->pd.DataFrame:
    '''
    Get a dataframe containing the fraction of cells_with specific one intID in one cluster
    '''
    fraction_cells_with_intID_in_cluster = []
    for i in range(df.shape[0]):
        fraction_cells_with_intID_in_cluster.append(
            df.loc[i, "num_cells_per_intID_in_cluster"] / df.loc[i, "num_cells_in_cluster"])
    df["fraction_cells_with_intID_in_cluster"] = fraction_cells_with_intID_in_cluster
    return df


def FilterStep1(df:pd.DataFrame,num_cell_in_cluster_cutoff:int,strict_cutoff_cells:float,strict_cutoff_fraction:float,fraction_cells_with_intID_in_cluster_cutoff:float) ->pd.DataFrame:
    '''
    Set threshold to filter the clusters with low quality or overwhelming noise
    '''
    num_cell_in_cluster_cutoff = num_cell_in_cluster_cutoff  # don't examine clusters with less than 5 cells
    strict_cutoff_cells = strict_cutoff_cells  # cutoff to define smaller clusters with more noise
    strict_cutoff_fraction = strict_cutoff_fraction  # for noiser, smaller clusters use a more strict "fraction_cells_with_intID_in_cluster_cutoff"
    fraction_cells_with_intID_in_cluster_cutoff = fraction_cells_with_intID_in_cluster_cutoff  # don't use intID to define a cluster if cluster has less than this fraction of cells with intID
    tmpdfToDivide_filtered = df[df.num_cells_in_cluster >= num_cell_in_cluster_cutoff] # Filter out clusters with a cell number of less than 5
    # For the clusters with a cell number of less than 20, the noise effect is more significant. For this small cluster, filter out
    # intIDs with a corresponding cell ratio of less than 0.35 (stricter threshold)
    tmpdfToDivide_filtered0 = tmpdfToDivide_filtered[
        (tmpdfToDivide_filtered.num_cells_in_cluster <= strict_cutoff_cells) & (
                    tmpdfToDivide_filtered.fraction_cells_with_intID_in_cluster >= strict_cutoff_fraction)]
    # For the clusters with a cell number of bigger than 20, filter out
    # intIDs with a corresponding cell ratio of less than 0.2 (stricter threshold)
    tmpdfToDivide_filtered1 = tmpdfToDivide_filtered[
        (tmpdfToDivide_filtered.num_cells_in_cluster >= strict_cutoff_cells) & (
                    tmpdfToDivide_filtered.fraction_cells_with_intID_in_cluster >= fraction_cells_with_intID_in_cluster_cutoff)]
    tmpdf_filtered_step1 = pd.concat([tmpdfToDivide_filtered0, tmpdfToDivide_filtered1])
    return tmpdf_filtered_step1



def Get_BarcodeSet(df:pd.DataFrame) ->pd.DataFrame:
    '''
    Get the static barcodes combination for each cluster
    '''
    cloneIDs_withdup=list(set(df.cloneID))
    intID_collections=[]
    for i in range(len(cloneIDs_withdup)):
        cloneID=cloneIDs_withdup[i]
        tmp=list(df[df.cloneID==cloneID].intID)
        tmp.sort()
        intID_collections.append(tmp)
    dic={}
    for i in range(len(cloneIDs_withdup)):
        dic[cloneIDs_withdup[i]]=intID_collections[i]
    cluster_int_df=pd.DataFrame(pd.Series(dic),columns=['intIDs']).reset_index().rename(columns={'index':'cloneID'})
    cluster_int_df=cluster_int_df.drop_duplicates("intIDs",keep="first")
    return cluster_int_df


def FilterStep2(df1:pd.DataFrame,df2:pd.DataFrame) ->pd.DataFrame:
    '''
    Remove the repeat clusters
    '''
    deduplicated_cloneIDs = list(df1.cloneID)
    deduplicated_index = []
    cloneIDlist = list(df2.loc[:, "cloneID"])
    for i in range(df2.shape[0]):
        if cloneIDlist[i] in deduplicated_cloneIDs:
            deduplicated_index.append(i)
    tmpdf_filtered_step2 = df2.iloc[deduplicated_index, :]
    return tmpdf_filtered_step2


def FilterStep3(df:pd.DataFrame) ->pd.DataFrame:
    '''
    To remove the incompletely same cluster,remove all clusters whose static barcodes combination is one subset of certain big clusters
    '''
    cluster_int_df=df[df.cloneID.isin(list(df.cloneID))]
    incomplete_repeat_clusters = ["Set"]*cluster_int_df.shape[0]
    for i in range(cluster_int_df.shape[0]):
        single_cluster_booleans=[""]*(cluster_int_df.shape[0]-1)
        temp=cluster_int_df[~cluster_int_df["cloneID"].isin([list(cluster_int_df.cloneID)[i]])]
        for j in range(cluster_int_df.shape[0]-1):
        #Check whether all static barcodes of cluster i are subset of static barcodes of cluster j
        #single_cluster_booleans means whether cluster i are subset of one cluster j
            count=0
            k=0
            while k<len(cluster_int_df.iloc[i,1]):
                if cluster_int_df.iloc[i,1][k] in temp.iloc[j,1]:
                    count=count+1
                k=k+1
            if count==len(cluster_int_df.iloc[i,1]):
                incomplete_repeat_clusters[i] = 'Subset' # Annotate subcluster as "subset"

    cluster_int_df["is_subset"]=incomplete_repeat_clusters
    cluster_int_df=cluster_int_df[cluster_int_df.is_subset!="Subset"] # Remove all subcusters
    #return the cluster * intID dataframe
    return cluster_int_df


def AssignCloneID(df1:pd.DataFrame,df2:pd.DataFrame,df3:pd.DataFrame) -> pd.DataFrame:
    '''
    #Identify the singlet/inter-clonal doublet and match/unmatch cells
    '''
    precluster_reads=df1
    intID_cluster=df2
    addcluster=df3
    unclustered_cells = precluster_reads[~precluster_reads.cellID.isin(list(addcluster.cellID))]
    #Add related information to the unclustered dataframe
    unclustered_cells0=unclustered_cells
    unclustered_cells=unclustered_cells0.copy()
    #How to assign cloneID/cluseterID to one cell?
    #1. Previous method
    #One intID may correspond to several clusters.
    #To solve this problem, the original paper merge the clusters that contain the repeated intIDs.
    #Thus, the intIDs in different clusters are completely different.
    #There is no intersection of any intIDs between different clusters.
    #Only when the intIDs of one cell are all belong to one cluster, the cell is defined as "singlet".

    #2. Method here
    #In our opinion, the combination of static intIDs are important for clone defintion.
    #If the intID combination between two clones are different, these two clones are different.
    # It's not necessary to ensure the intIDs are completely different from each others among different clones.
    #The method in the paper merge small clones into a big clones which contains many intIDs.
    # The considerable number of intIDs in this big clones may make this clone contain too many cells.
    #Thus, for the cluster assignment issue, our method only check whether the all intIDs of one cell belong to one specific clone.
    #If all intIDs of this cell belong to one cluster, this cell will be assinged as this cloneID.
    #Since we delete all the clones that contains subset intID combination, each cell belongs to at most one cluster.
    #The cells whose intIDs do not belong to any clusters will be defined as "unmatched".


    cells = list(unclustered_cells.loc[:,'cellID'].unique())
    all_barcodes = list(intID_cluster.loc[:,'intID'].unique())
    #Get all unique barcode intIDs in clustered and unclustered cells
    all_barcodes = all_barcodes+list(unclustered_cells.loc[:,'intID'].unique())
    all_barcodes = list(set(all_barcodes))
    barcodes_cells = all_barcodes
    #Get all unique barcode clusterIDs/cloneIDs in clustered cells
    clusters = list(intID_cluster.loc[:,'cloneID'].unique())
    barcodes_clusters = all_barcodes
    #Build the empty unique barcode * cell matrix and unique barcode * cloneID matrix
    cell_barcode = pd.DataFrame(index=barcodes_cells, columns=cells)
    cell_barcode.fillna(0, inplace=True)
    cluster_barcode = pd.DataFrame(index=barcodes_clusters, columns=clusters)
    cluster_barcode.fillna(0, inplace=True)
    #Genetare the corresponding intID-cell information in the matrix
    for i in range(unclustered_cells.shape[0]):
        cell_barcode.loc[unclustered_cells.iloc[i].loc['intID'],unclustered_cells.iloc[i].loc['cellID']] = 1
    #Genetare the corresponding intID-cloneID information in the matrix
    for i in range(intID_cluster.shape[0]):
        cluster_barcode.loc[intID_cluster.iloc[i].loc['intID'],intID_cluster.iloc[i].loc['cloneID']] = 1
    cell_cluster = pd.DataFrame(columns=['cellID','cloneID'])
    cell_infor=[]
    match_infor=[]
    singlet_infor=[]
    status_infor=[]
    clone_infor=[]
    #Screen for all cells in the unclustered cells
    for i in range(len(cells)):
        belongto = []
        #Get the all intIDs for cells i
        codeofcell = cell_barcode.loc[:,cells[i]]
        is_match="unmatch"
        is_singlet="inter-clonal doublets"
        status=""
        for x in range(len(clusters)):
            #Screen for all clusters in clustered cells. Here get the intIDs for cluster x
            codeofcluster = cluster_barcode.loc[:,clusters[x]]
            keep = True
            #Screen all intIDs for cell i
            for y in codeofcell.index:
                #Check if there is intID of this cell correspond to one cluster
                if codeofcell[y] ==1 and codeofcluster[y]==1:
                    is_match="match"
                #Screen for all intID in cell i and cluster x.
                #If there is one intID for cell i that does not belong to cluster x, cell i will not be assigned to cluster x
                #Only when all intIDs of one cell are subset of intIDs of one cluster, this cell will be assigned to this cluster
                if codeofcell[y] ==1 and codeofcluster[y]==0:
                    keep = False
            if keep:
                #cell i belong to cluster x
                belongto.append(clusters[x])
        #The belongto store the cloneIDs this cell corresponds.
        if len(belongto) == 1:
            is_singlet="singlet"
        if is_match=="match":
            if is_singlet=="singlet":
                status="singlet"
            if is_singlet=="inter-clonal doublets":
                status="inter-clonal doublets"
        else:
            status="unmatch"
        if len(belongto)==0:
            belongto.append("No singlet")
        #is_match,is_singlet,assigned ID and status are collected
        cell_infor.append(cells[i])
        match_infor.append(is_match)
        singlet_infor.append(is_singlet)
        status_infor.append(status)
        clone_infor.append(belongto[0])
    #Add all information to the matrix of the unclustered cells
    result=pd.DataFrame()
    result["cellID"]=cell_infor
    result["is_match"]=match_infor
    result["is_singlet"]=singlet_infor
    result["status"]=status_infor
    result["cloneID"]=clone_infor


    if unclustered_cells.empty and result.empty:
        All_withcloneID = precluster_reads
        All_withcloneID["is_match"]=["match"]*All_withcloneID.shape[0]
        All_withcloneID["is_singlet"]=["singlet"]*All_withcloneID.shape[0]
        All_withcloneID["status"]=["singlet"]*All_withcloneID.shape[0]
        All_withcloneID=precluster_reads.merge(addcluster,on="cellID")



    else:
        unclustered_cells_cloneID=unclustered_cells.merge(result,on="cellID",how="left")
        #Add all information to the matrix of the clustered cells
        clustered_cells = precluster_reads[precluster_reads.cellID.isin(list(addcluster.cellID))]
        clustered_cells["is_match"]=["match"]*clustered_cells.shape[0]
        clustered_cells["is_singlet"]=["singlet"]*clustered_cells.shape[0]
        clustered_cells["status"]=["singlet"]*clustered_cells.shape[0]
        clustered_cells=clustered_cells.merge(addcluster,on="cellID")
        All_withcloneID=pd.concat([clustered_cells,unclustered_cells_cloneID])
    return All_withcloneID


def Clones_in_PT(df:pd.DataFrame,primary_position:str) ->pd.DataFrame:
    '''
    Get the clones that exist in the primary tumor and calculate the
    '''

    if primary_position in list(set(df.Sample_posi)):
        print("Get the cloneIDs in primary tumor")
        forsingletstep3=df[df.Sample_posi==primary_position]
    else:
            print("Primary position does not exist and Clone_in_PT step is skipped")
            forsingletstep3=df

    return forsingletstep3


def FilteredClone(precluster_reads:pd.DataFrame,addcluster:pd.DataFrame,primary_position:str) ->pd.DataFrame:
    ##########################################################
    print("Prepare for filter step")
    precluster_reads1 = precluster_reads
    ToFilter = precluster_reads[precluster_reads.cellID.isin(list(addcluster.cellID))]
    ToFilter = ToFilter.merge(addcluster, on="cellID")
    UniqueToFilter = ToFilter.drop_duplicates(['cellID', "intID"], keep='first')
    # Get the number of cells in each cluster
    Clusters_count = Get_num_of_cells_in_each_cluster(UniqueToFilter)
    # Get the number of cells that each intID corresponds in each cluster
    tmpdf = Get_num_cells_per_intID_in_cluster(UniqueToFilter)
    tmpdfToDivide = tmpdf.merge(Clusters_count, on="cloneID")
    # Get the fraction of cells that one intID correspond in one cluster
    tmpdfToDivide = Get_fraction_cells_with_intID_in_cluster(tmpdfToDivide)
    # Filter the cluster accroding to the number of cells in this cluster
    # Filter the intID that used to define the cluster according to the fraction of cells that contain this intID
    tmpdf_filtered_step1 = FilterStep1(tmpdfToDivide, 5, 20, 0.35, 0.2)
    print("Filter the clusters and intIDs with low quality")
    # Get the intID combination in each cluster
    cluster_int_df = Get_BarcodeSet(tmpdf_filtered_step1)
    print("Get the intID combination for each clone")
    # Remove the repeat clusters(same intID combination)
    tmpdf_filtered_step2 = FilterStep2(cluster_int_df, tmpdf_filtered_step1)
    # Remove all subclusters(the intID combination of subcluster is the subset of a bigger cluster)
    # Get a dataframe that use intIDs to define a cluster or clone
    cluster_int_df = FilterStep3(cluster_int_df)
    # Get the all intID+cellID information for each cells
    cell_intID_df = precluster_reads1
    # Get the all cloneID information related each intID
    cluster_intID_df = tmpdf_filtered_step2.loc[:,["intID","cloneID"]]
    # Assign cloneID to all cells
    print("Assign a cloneID to each cell")
    forsingletstep1=AssignCloneID(precluster_reads1,cluster_intID_df,addcluster)
    print("Extract the singlet for downstream analysis")
    forsingletstep2 = forsingletstep1[forsingletstep1.status == "singlet"]
    forsingletstep2["Sample_posi"]=forsingletstep2["sample"]
    forsingletstep3 = Clones_in_PT(forsingletstep2,primary_position)
    cloneIDs_in_PT = list(set(forsingletstep3.cloneID))
    # Preserve the cells from the clones existing in the primary tumors
    Final = forsingletstep2[forsingletstep2.cloneID.isin(cloneIDs_in_PT)]
    print("Filter step finish!")
    return Final

