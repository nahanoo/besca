#this module contains functions for cell type annotation based on signatures in python using scanpy
#import using the python version 1.3.2 at least is a prerequisit! needs to be checked in all functions

from scanpy.api import AnnData
from pandas import DataFrame, read_csv
from scipy.stats import mannwhitneyu
from numpy import log, where


def getset(df, signame, Cell, thresholdTable):
    """ Handles missing signatures aux function for make_anno
    Based on a dataframe of p-values, a signature name and a cutoff check if sign is present
    parameters
    ----------
    df: panda.DataFrame
      a dataframe of p-values per signature and cell cluster
    signame: str
      signature prefix
    Cell: str
      celltype to observe. Should match key in threshold Table and one suffix in df columns
    threshold: numpy.float64
      cutoff used for cluster attributions
    returns
    -------
    set
        Set of clusters passing the threshold value
    """
    signame_complete = signame+'_' + Cell
    if(bool(df.index.isin([signame_complete]).any())):
        if(Cell in thresholdTable.keys()):
            return(set(df.columns[df.loc[signame_complete, :] > thresholdTable[Cell]['0']]))
        else:
            # Should we add verbosity
            return(set(df.columns[df.loc[signame_complete, :] > thresholdTable['default_value']['0']]))
    else:
        return(set())


def gene_fr(f, Tcells):
    """ Returns the fraction of cells expressing a gene in distinct clusters
    Takes as an input a dataframe with fractions per clusters and a gene name
    Returns the fraction genes per cluster
    parameters
    ----------
    f: panda.DataFrame
      a dataframe of fraction genes per cell as outputted by besca
    Tcells: 'str'
      gene name

    returns
    -------
    panda.Series
        a Series of fractions per cluster and gene
    """
    k = f.loc[f['Description'].isin( Tcells),].copy()
    k = k.iloc[:,0:len(k.columns)]
    return(k.mean(axis=0, skipna=True))


def score_mw(f, mymarkers):
    """ Score Clusters based on a set of immune signatures to generate a df of pvals
    Takes as an input a dataframe with fractions per clusters and a dictionary of signatures
    Performs a Mann-Whitney test per each column and Signature, returns -10logpValues
    parameters
    ----------
    f: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    mymarkers: Dictionary
      a dictionary of signatures

    returns
    -------
    panda.DataFrame
        a dataframe of -10logpValues per cluster and signature
    """
    ids = set(f.iloc[:,2:len(f.columns)])
    mypFrame = DataFrame(index=mymarkers.keys(), columns=ids)
    for key, value in mymarkers.items():
        for i in ids:
            mypFrame.loc[key,i] = -10*log(mannwhitneyu(x=f.loc[f['Description'].isin(value),:][i], 
                    y=f[i], alternative='greater').pvalue)
    return(mypFrame)




def add_anno(adata,cNames, cluster='louvain'):
    """ Adds annotation generated with make_anno to a AnnData object
    Takes as input the AnnData object to which annotation will be appended and the large annotation
    Generates 5 new adata.obs columns with annotation at distinct levels
    parameters
    ----------
    adata: AnnData
      AnnData object that is to be annotated
    cNames: panda.DataFrame
      a list of cluster names generated with make_anno
    cluster: string
      initial cluster used for annotation (column in adata.obs)

    returns
    -------
    AnnData
        the annotated adata object
    """
    adata.obs['clusterID'] = adata.obs[cluster]
    adata.rename_categories('clusterID', cNames)
    
    ### split according to T, B, monocytes, tumor cells
    cellnames = []
    cellgroups = []
    scellgroups = []
    sscellgroups = []
    for i in adata.obs['clusterID']:
        cellnames.append(i.split(".")[0]+i.split(".")[2]+i.split(".")[3])
        cellgroups.append(i.split(".")[2])
        scellgroups.append(i.split(".")[2]+i.split(".")[3])
        sscellgroups.append(i.split(".")[2]+i.split(".")[3]+i.split(".")[4])
    adata.obs['cell_names'] = cellnames
    adata.obs['cell_group'] = cellgroups
    adata.obs['scell_group'] = scellgroups
    adata.obs['sscell_group'] = sscellgroups
    return(adata)


def make_anno(mypFrame, signame, f, thresholdDict,
              CD45threshold=0.3, species='human'):
    """ Annotate Immune cells and some other cell types (melanoma, endothelial)
    Based on a dataframe of -log10pvals, a cutoff and a signature set generate cell annotation
    Hierarchical model of Immune cell annotation.
    It expects a specific set of signatures with a common prefix (signame) and specified suffixes indicating
    the cell type.
    parameters
    ----------
    mypFrame: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    f: panda.DataFrame
      a dataframe with fraction genes expressed per cluster
    signame: 'str'
      signature base name; should be prefix of multiples columns in mypFrame
    thresholdDict: 'dict'
      dictionary containing the thresholds. expected: key: cell name 
      value: threshold.   See function generate_annot_threshold and get_threshold
    CD45threshold: 'float'
      threshold for CD45 percent. Cluster with higher percent of cell expressing of CD45 (PTPRC/Ptprc) 
      will be assigned as possible immune cell.

    returns
    -------
    list of str
        String of cluster labels
    """
    if(not 'default_value' in thresholdDict.keys()):
        maxI = max(df['0'])
        thresholdDict['default_value'] = {'0': maxI}
    ## over CD45threshold % of cells have CD45
    aCD45P = set(where(gene_fr(f, ['PTPRC']) >= CD45threshold)[0])
    if species == 'mouse':
        aCD45P = set(where(gene_fr(f, ['Ptprc']) >= CD45threshold)[0])
    aBcells = getset(mypFrame, signame, 'Bcells', thresholdDict)
    aPlasma = getset(mypFrame, signame,'Plasma', thresholdDict)
    aTcd8 = getset(mypFrame, signame, 'Tcd8', thresholdDict)
    aTcd4 = getset(mypFrame, signame, 'Tcd4', thresholdDict)
    aTcgd = getset(mypFrame, signame, 'Tcgd', thresholdDict)
    aTreg = getset(mypFrame, signame, 'Treg', thresholdDict)
    aTNK = getset(mypFrame, signame, 'TNK', thresholdDict)
    aTilCM = getset(mypFrame, signame, 'TilCM', thresholdDict)
    aT4CM = getset(mypFrame, signame, 'T4CM', thresholdDict)
    aTcytox = getset(mypFrame, signame, 'Tcytox', thresholdDict)
    aTEM = getset(mypFrame, signame, 'TEM', thresholdDict)
    aTtexh = getset(mypFrame, signame, 'Ttexh', thresholdDict)
    aTpexh = getset(mypFrame, signame, 'Tpexh', thresholdDict)
    aNKT = getset(mypFrame, signame, 'NKT', thresholdDict)
    aTcells = getset(mypFrame, signame, 'Tcells' ,thresholdDict)
    aNKcells = getset(mypFrame, signame, 'NKcells', thresholdDict)
    aNKnai = getset(mypFrame, signame, 'NKnai', thresholdDict)
    aNKcyt = getset(mypFrame, signame, 'NKcyt', thresholdDict)
    aMyelo = getset(mypFrame, signame, 'Myelo', thresholdDict)
    aNai = getset(mypFrame, signame, 'Naive', thresholdDict)
    aAct = getset(mypFrame, signame, 'Activation', thresholdDict)
    aMem = getset(mypFrame, signame, 'Memory', thresholdDict)
    aEff = getset(mypFrame, signame, 'Eff', thresholdDict)
    aNonEff = getset(mypFrame, signame, 'NonEff',thresholdDict)
    aCytotox = getset(mypFrame, signame, 'Cytotox',thresholdDict)
    aDCs = getset(mypFrame, signame, 'aDCs', thresholdDict)
    moDC = getset(mypFrame, signame, 'moDC', thresholdDict)
    pDCs = getset(mypFrame, signame, 'pDCs', thresholdDict)
    acDC1 = getset(mypFrame, signame, 'cDC1', thresholdDict)    
    acDC2 = getset(mypFrame, signame, 'cDC2', thresholdDict)   
    aMono = getset(mypFrame, signame, 'Monocytes', thresholdDict)
    aTAM = getset(mypFrame, signame, 'TAM', thresholdDict)
    aTMO = getset(mypFrame, signame, 'TMO', thresholdDict)
    aMo14 = getset(mypFrame, signame, 'Mo14', thresholdDict)
    aMo16 = getset(mypFrame, signame, 'Mo16', thresholdDict)
    aTAMCx = getset(mypFrame, signame, 'TAMCx', thresholdDict)
    aTMid = getset(mypFrame, signame, 'TMid', thresholdDict)
    aMoMa = getset(mypFrame, signame, 'MoMa', thresholdDict)
    aMyeloSubtype = getset(mypFrame, signame, 'MyeloSubtype', thresholdDict)
    aGranu = getset(mypFrame, signame, 'Granulo', thresholdDict)
    aNeutro = getset(mypFrame, signame, 'Neutrophil', thresholdDict)
    aMacrophage = getset(mypFrame, signame, 'Macrophage', thresholdDict)
    aCC = getset(mypFrame, signame, 'Cellcycle', thresholdDict)
    aCh = getset(mypFrame, signame, 'Checkpoint', thresholdDict)
    aIfng = getset(mypFrame, signame, 'Ifng', thresholdDict)
    aEndo = getset(mypFrame, signame, 'Endo', thresholdDict)
    aCafs = getset(mypFrame, signame, 'Cafs', thresholdDict)
    aMelMelan = getset(mypFrame, signame, 'MelMelan', thresholdDict)
    aMega = getset(mypFrame, signame, 'Megakaryocytes', thresholdDict)
   #### Part 2 Combine
    CD45P = set(aCD45P).union(set(aBcells).union(aNKcells).union(aTNK).union(aTcells).union(aTreg).union(aMyelo).union(pDCs))
    Tc = set(CD45P).intersection(set(aTcells).union(aTreg)-set(aBcells)-set(aMyelo))
    Nk = set(CD45P).intersection(set(aNKcells).union(set(aCytotox))-set(aTcells))
    Pl = set(aPlasma).intersection(set(CD45P))-Tc-Nk
    Bc = set(aBcells).union(Pl)
    Tr = set(Tc).intersection(aTreg)
    Dc = set(aDCs).union(acDC1).union(acDC2).union(moDC).union(pDCs)-Tc-Nk-Bc-Pl
    My = set(CD45P).intersection(set(aMoMa).union(set(aTAM).union(set(aTMO)).union(set(aTMid)))).union(Dc)-Tc-set(aBcells)-Nk
    naiTc = set(Tc).intersection(set(aNai))
    ccTc = set(Tc).intersection(set(aCC))
    chTc = set(Tc).intersection(set(aCh))
    effTc = set(Tc).intersection(set(aCytotox)-Tr)
    actTc = set(Tc).intersection(set(aAct)-Tr)
    memTc = set(Tc).intersection(set(aMem))-actTc-Tr-naiTc
    exTc = set(Tc).intersection(chTc.union(aNonEff))-naiTc
    texTc = set(Tc).intersection(aTtexh)-naiTc
    pexTc = set(Tc).intersection(aTpexh.intersection(chTc))-naiTc
    cytoTc = set(Tc).intersection(set(aTcytox)-Tr)-naiTc
    NKTc = set(Tc).intersection(set(aNKT)-Tr)
    TilCM = set(Tc).intersection(set(aTilCM)-Tr) #IL7R cluster
    T4CM = set(Tc).intersection(set(aT4CM)-Tr)-TilCM-cytoTc
    TEM = set(Tc).intersection(set(aTEM).intersection(aTcytox)-Tr).union(aEff)-TilCM-T4CM-cytoTc
    TcCD8 = set(aTcd8).intersection(Tc)
    Tr = Tr-TcCD8
    TcCD4 = set(aTcd4).intersection(Tc)-Tr-set(TcCD8)
    Tgd = set(Tc).intersection(aTcgd)
    Nknai = aNKnai.intersection(Nk)
    Nkcyt = aNKcyt.intersection(Nk)-Nknai
    Mo14 = aMo14.intersection(My)-acDC1-acDC2-aDCs-Tc-Bc-Nk
    Mo16 = aMo16.intersection(My)-acDC1-acDC2-aDCs-Tc-Bc-Nk
    Mo14 = Mo14-Mo16
    pDCs = pDCs-Bc
    #### Part 3 annotate
    ids = list(set(f.iloc[:, 2:len(f.columns)]))
    nclust = len(ids)
    ###### Summarize the annotations:
    cNames = []
    CD45P = [str(i) for i in CD45P]
    for i in range(0,nclust):
        ii = str(i)
        newname = "C"+ii+['.nCD45', '.CD45'][(ii in CD45P)*1]
        if (ii not in CD45P):
            lev2 = ['no', '.Endo'][(ii in aEndo)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.Mel'][(ii in aMelMelan)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.Caf'][(ii in aCafs)*1]
            lev3='.'
            lev4='.'
            
        if (ii in CD45P):
            lev2 = ['no', '.Tc'][(ii in Tc)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.NK'][(ii in Nk)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.Bc'][(ii in Bc)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.My'][(ii in My)*1]
            if (lev2 == 'no'):
                lev2='.'
                lev3='.'
                lev4='.'
    
        if (lev2 == '.Tc'):
            lev3 = ['no', '.8'][(ii in TcCD8)*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Tr'][(ii in Tr)*1]
            if (lev3 == 'no'):
                lev3 = ['.4', '.4'][(ii in TcCD4)*1]
            lev4 = ['no', '.Tr'][(ii in Tr)*1]          
            if (lev4 == 'no'):
                lev4 = ['no', '.texh'][(ii in texTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.NKT'][(ii in NKTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Tgd'][(ii in Tgd)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Cytox'][(ii in cytoTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.CMil7'][(ii in TilCM)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.CM4'][(ii in T4CM)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.CC'][(ii in ccTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.texh'][(ii in pexTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.EM'][(ii in TEM)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Nai'][(ii in naiTc)*1]  
            if (lev4 == 'no'):
                lev4 = ['no', '.Exh'][(ii in chTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Eff'][(ii in effTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Act'][(ii in actTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Mem'][(ii in memTc)*1]
            if (lev4 == 'no'):
                lev4 = '.Mid'
        if (lev2 == '.Bc'):
            lev3 = ['.Bc', '.Pl'][(ii in Pl)*1]
            lev4 = '.'
        if (lev2 == '.NK'):
            lev3 = ['no', '.Nai'][(ii in Nknai)*1]
            if (lev3 == 'no'):
                lev3 = '.Cyt'
            lev4 = '.'               
        if (lev2 == '.My'):
            lev3 = ['no', '.pDC'][(ii in pDCs.intersection(My))*1]   
            if (lev3 == 'no'):
                lev3 = ['no', '.aDCs'][(ii in aDCs.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.cDC1'][(ii in acDC1.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.cDC2'][(ii in acDC2.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.aTc'][(ii in aTcells.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Mo14'][(ii in Mo14.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Mo16'][(ii in Mo16.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.TAMCx'][(ii in aTAMCx.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.TMid'][(ii in aTMid.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.aTAM'][(ii in aTAM.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.moDC'][(ii in moDC.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.TMO'][(ii in aTMO.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.MoMa'][(ii in aMyeloSubtype.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Gra'][(ii in aGranu.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Mky'][(ii in aMega.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Neu'][(ii in aNeutro.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Tc'][(ii in aTcells)*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.My'][(ii in aMyelo)*1]
            if (lev3 == 'no'):
                lev3 = '.'
            lev4 = '.'
        newname = newname+lev2+lev3+lev4+['.nNa', '.Na'][(ii in aNai)*1]
        newname += ['.nCC', '.CC'][(ii in aCC)*1]+['.nCy', '.Cy'][(ii in aCytotox)*1]
        newname += ['.nCh', '.Ch'][(ii in aCh)*1]+['.nAc', '.Ac'][(ii in aAct)*1]
        cNames.append(newname)
    return(cNames)




