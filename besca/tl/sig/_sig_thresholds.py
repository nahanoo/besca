import pandas as pd

def generate_annot_threshold(myc, filename='sig_threshold.csv'):
    """ Generate a config files for autoannotation of the clusters.
    This file will be used afterwards by the sig modules.
    Modifying the thresholds is then possible on the text file
    to take into account datasets specificity.
    parameters
    ----------
    myc: 'numpy.float64'
      threshold used for cluster attribution
    fileName: 'str'
      where this will be saved.
    returns
    -------
     None
    Example
    -------
    >>> import besca as bc 
    >>> bc.tl.sig.generate_annot_threshold(3)

    """
    cellName = pd.Series(['Bcells', 'Plasma', 'Tcd8', 'Tcd4',
                'Tcgd', 'Treg', 'TNK', 'TilCM', 'T4CM',
                'Tcytox', 'TEM', 'Ttexh', 'Tpexh', 
                'NKT', 'Tcells', 'NKcells', 'NKnai',
                'NKcyt', 'Myelo', 'Naive', 'Activation',
                'Memory', 'Eff', 'nonEff', 'Cytotox',
                'aDCs', 'moDC', 'pDCs', 'cDC1', 
                'cDC2', 'Monocytes', 'TAM', 'TMO',
                'Mo14', 'Mo16', 'TAMCx', 'TMid',
                'MyeloSubtype', 'Granulo',  'Neutrophil', 'Macrophage',
                'Cellcycle'  , 'Checkpoint', 'Ifng',
                'Endo', 'Cafs',  'MelMelan', 'Megakaryocytes'])

    threshold = pd.Series([1, 1, 1, 1/2, 2/3,
                 1, 2, 2, 3/2,
                 2, 2, 1, 1, 3/2, 4/3, 2, 1,#Nknai
                 2, 1, 1, 1, 1, 1, 1, 2,#Cytotox
                 1, 3/2, 2, 1, 2, 1, 3/2, 1,#TMO
                 1, 1, 1, 2,# Moma
                 1, 1, 1, 1,# Macrophage
                 4/3, 1, 1,#Ifng
                 1, 1, 1, 1])# Megakaryocytes
    threshold = pd.Series([x * myc for x in threshold], index = cellName)
    threshold.to_csv(filename)


def get_threshold(filethreshold='sig_threshold.csv'):
    """ Read a threshold file,csv, and return a dictionary.
    Initial column contains keys an second a value.
    If not present default_value key will be added 
    ----------
    myc: 'numpy.float64'
      threshold used for cluster attribution
    fileName: 'str'
      where this will be saved.
    returns
    -------
    a dict. 
        Key: first columns of the csv file
        Value: Dictionnary : {'0' : value}
            where the value is the second column of the csv file
    Example
    -------
    >>> import besca as bc 
    >>> bc.tl.sig.generate_annot_threshold(3, filename='test.csv')
    >>> thresholdDict = bc.tl.sig.get_threshold('test.csv')
    """
    df = pd.read_csv(filethreshold, index_col=0)
    myc_dict = df.to_dict(orient='index')
    if(not 'default_value' in myc_dict.keys()):
        maxI = max(df['0'])
        myc_dict['default_value'] = {'0': maxI}
    return(myc_dict)



