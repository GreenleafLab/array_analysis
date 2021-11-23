import os
import sys
import numpy as np
import pandas as pd
import pickle
import itertools
import datetime

def returnFigDirectory():
    return 'figs_%s'%str(datetime.date.today())

def stripExtension(filename):
    end_idx = min([idx for idx in [filename.find('.pkl'), filename.find('.gz'), len(filename)] if idx >= 0])
    return os.path.splitext(filename[:end_idx])[0]

def loadFile(filename):
    """ Find extension and return loaded file. """
    if filename is None:
        print("Error: No filename given!")
        sys.exit()
    ext = os.path.splitext(filename)[-1]
    
    # check compression status
    if ext == '.gz':
        compression = 'gzip'
        ext = os.path.splitext(filename[:filename.find('.gz')])[-1]
    else:
        compression = None

    if ext == '.json.zip':
        return pd.read_json(filename)
    
    elif ext == '.pkl':
        return pd.read_pickle(filename)

    elif ext == '.CPseq':
        return _loadCPseq(filename, compression=compression)
    
    elif ext == '.BCseq' or ext == '.libChar':
        return _loadUnindexedTable(filename, compression=compression)
    
    elif ext == '.CPfluor':
        return _loadCPFluorFile(filename, compression=compression)

    elif ext == '.CPvariant' or ext == '.CPseries' or ext == '.CPfitted' or ext == '.fitParameters' or ext == '.CPannot':
        return _loadIndexedTable(filename)
   
    elif ext == '.txt':
        return _loadTextFile(filename)
    
    elif ext == '.p':
        return _loadPickle(filename)
    
    elif ext == '.times':
        return _loadTextFile(filename)

    else:
        print('Extension %s not recognized. No file loaded.'%ext)
        
def saveFile(filename, data, **kwargs):
    """Save data to a file according to extension."""
    if filename is None:
        print("Error: No filename given!")
        sys.exit()
    ext = os.path.splitext(filename)[-1]
    if ext == '.pkl':
        data.to_pickle(filename)
    elif ext == '.CPvariant' or ext == '.dat':
        data.to_csv(filename, sep='\t')
    elif ext == '.p':
        pickle.dump(data, open( filename,  "wb" ))
    elif ext == '.csv':
        data.to_csv(filename, **kwargs)

    else:
        print('Extension %s not recognized. No file saved.'%ext)
        
    
def _loadCPseq(filename, **kwargs):
    """ Return CPseq file. """
    return pd.read_table(filename, header=None,
                         names=['clusterID', 'filter',
                                'read1_seq', 'read1_quality',
                                'read2_seq', 'read2_quality',
                                'index1_seq','index1_quality',
                                'index2_seq', 'index2_quality'],
                         index_col=0, **kwargs)

def _loadUnindexedTable(filename, **kwargs):
    """ Return tab separated table """
    return pd.read_table(filename, **kwargs)

def _loadIndexedTable(filename, **kwargs):
    """ Return tab separated table with an index col """
    return pd.read_table(filename, index_col=0, **kwargs)

def _loadCPFluorFile(filename, **kwargs):
    a = pd.read_csv(filename,  usecols=list(range(7, 12)), sep=':', header=None, names=['success', 'amplitude', 'sigma', 'fit_X', 'fit_Y'], **kwargs )
    b = pd.read_csv(filename,  usecols=list(range(7)), sep=':', header=None,  dtype=str, **kwargs)
    a.index = (b.loc[:, 0] + ':' + b.loc[:, 1] + ':' + b.loc[:, 2] + ':' +
               b.loc[:, 3] + ':' + b.loc[:, 4] + ':' + b.loc[:, 5] + ':' + b.loc[:, 6])
    a.index.name = 'clusterID'
    return a

def _loadTextFile(filename):
    try:
        a = np.loadtxt(filename)
    except ValueError:
        a  = np.loadtxt(filename, dtype=str)
    return a

def _loadPickle(filename):
    return pickle.load( open( filename, "rb" ) )


