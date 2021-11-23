import numpy as np
import scipy.stats as st
import pandas as pd


def reverseComplement(seq, rna=None):
    """
    This is a function that gives the reverse complement of any sequence
    """
    if rna is None: rna = False
    for base in seq:
        if base not in 'NATCGUatcgu':
            print "TypeError: NOT a DNA sequence\t%s"%seq
            return None
    seq_dict = {}
    seq_complement = 'TAGCAN'
    
    # if upper case
    for i, base in enumerate('ATCGUN'):
        seq_dict[base] = seq_complement[i]
    
    # if lower case
    for i, base in enumerate('atcgu'):
        seq_dict[base] = seq_complement[i].lower()
    if rna:
        return "".join([seq_dict[base] for base in reversed(seq)]).replace('T', 'U')
    else:
        return "".join([seq_dict[base] for base in reversed(seq)])

def rc(seq, rna=None):
    return reverseComplement(seq, rna=rna)

def poisson(k, lambda_):
    
    return np.power(lambda_,k)*np.exp(-lambda_)/scipy.misc.factorial(k)

def getFDR(score, nullscores):
    """
    use null distribution to get false discovery rate of scores
    """

    if score > np.median(nullscores):
        q = np.sum(nullscores > score)/float(len(nullscores))
    elif score < np.median(nullscores):
        q = np.sum(nullscores < score)/float(len(nullscores))
    else:
        q = 1
            
    return q

def getFDRs(scores, nullscores):
    """
    use null distribution to get false discovery rate of scores
    """
    q = np.zeros(len(scores))
    for i, score in enumerate(scores):
        q[i] = getFDR(score, nullscores)
            
    return q

def getFDR_onetailed(scores, nullscores):
    """
    use null distribution to get false discovery rate of scores
    """
    q = np.zeros(len(scores))
    for i, score in enumerate(scores):
        if i%1000==0: print 'iteration %d'%i
        q[i] = np.sum(np.abs(nullscores) > np.abs(score))/float(len(nullscores)) 
    return q

def getCDF(x):
    """
    returns the bins and cumulative frequencis of vector 'x'
    """
    xvalues = np.sort(x)
    yvalues=np.arange(len(xvalues))/float(len(xvalues))
    return xvalues, yvalues


def divideIntoBins(bin_by_me, bin_me, numBins=None, binEdges=None):
    """
    use variable 'bin_by_me' to divide into numBins number of bins.
    Find the values of bin_me that are in each of these bins.
    Also allows you to inout the binedges yourself rather than finding them by binning 'bin_by_me'
    """
    if numBins is None:
        numBins = 20
    if binEdges is None:
        binEdges = np.array(np.percentile(bin_by_me, np.linspace(0, 100, numBins+1).tolist()))
    else: numBins = len(binEdges)-1
    binCenters = (binEdges[:-1] + binEdges[1:])*0.5
    reads2binIndxRight = np.digitize(bin_by_me, bins=binEdges, right=True)   # for every peak, this is telling you the bin it belongs in
    
    binIndxRight = np.arange(numBins)+1 # these are the indexes of the bins (unique)
    
    vecs = ['']*numBins
    for binIndxLeft in range(numBins):
        
        indxSubsetInBin = reads2binIndxRight == binIndxRight[binIndxLeft]
        
        vecs[binIndxLeft] = bin_me[indxSubsetInBin]
        # get rid of nans and inf because those suck

    return vecs, binEdges

def is_outlier(points, thresh=None):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if thresh is None:
        thresh = 3.5
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def remove_outlier(points, thresh=None):
    return points[np.logical_not(is_outlier(points, thresh=thresh))]

def getCorrelation(vec1, vec2):
    index = np.all(np.isfinite(np.vstack([vec1, vec2]).astype(float)), axis=0)
    return st.pearsonr(vec1[index], vec2[index])[0]

def fillNAMat(submat):
    """For any NaN values in matrix, fill with the average value for that column."""
    submatnew = submat.copy()
    for col in submat:
        submatnew.loc[submat.loc[:, col].isnull(), col] = submat.loc[:, col].mean()
    return submatnew


def doPCA(submat, fillna=False):
    """Perform PCA on a matrix."""
    if fillna:
        submat = fillNAMat(submat)

    from fittinglibs import dummylib
    return dummylib.doPCA(submat)

def transform_data(data, loadings, scale_params=None):
    """Transform data in coordinates given by 'loadings'.

    If 'scale_params' is given, subtract off mean and divide by std dev of data first."""
    if scale_params is not None:
        failed = True
        if len(scale_params)==2:
            if len(scale_params[0])==data.shape[1] and len(scale_params[1]==data.shape[1]):
                data = ((data - scale_params[0])/scale_params[1]).copy()
                failed = False
        if failed:
            print 'check dimensions of scale_params'
    
    transformed = pd.DataFrame(np.dot(loadings, data.transpose())).transpose()
    transformed.index = data.index
    transformed.columns = ['pc_%d'%i for i in np.arange(transformed.shape[1])]
    return transformed

