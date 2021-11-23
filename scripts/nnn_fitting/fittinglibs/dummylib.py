import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.misc
import scipy.stats as st
import pandas as pd
from sklearn.decomposition import PCA as sklearnPCA


def doPCA(submat):
    """Perform PCA on a matrix."""
    sklearn_pca = sklearnPCA(n_components=None, whiten=False)
    sklearn_transf = sklearn_pca.fit_transform(submat)
    pca_projections = np.dot(np.linalg.inv(np.dot(sklearn_transf.T, sklearn_transf)),
                                 np.dot(sklearn_transf.T, submat))
    
    
    return (sklearn_pca, pd.DataFrame(sklearn_transf, index=submat.index, columns=['pc_%d'%i for i in np.arange(sklearn_transf.shape[1])]),
            pd.DataFrame(pca_projections, columns=submat.columns, index=['pc_%d'%i for i in np.arange(sklearn_transf.shape[1])]))