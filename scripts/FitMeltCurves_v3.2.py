import pandas as pd
import numpy as np
import sys
import scipy.optimize as opt
from tqdm import tqdm
import argparse
from pandarallel import pandarallel
from datetime import date
from glob import glob
import os

todaysdate = date.today().strftime("%d%b%Y")
pandarallel.initialize(progress_bar=True)
tqdm.pandas()

def clip_normalize(data):
    '''
    Used when normalizing signal.
    '''

    min_fluor = np.nanpercentile(data,1)
    max_fluor = np.nanpercentile(data,99)
    #max_fluor = np.nanmedian(data)+5*np.nanstd(data)
    return np.clip(data, min_fluor, max_fluor)

def percentile(n):
    '''
    Helper to return percentile in pandas aggregate call.
    '''
    def percentile_(x):
        return np.nanpercentile(x,n)
    percentile_.__name__='percentile_%s' % n
    return percentile_

def bootstrap_inds(n):
    return np.random.choice(n,n)

def resample_df(df):
    bs_inds = bootstrap_inds(len(df))
    new = df.iloc[bs_inds]
    new = new.reset_index()
    return new

def bootstrap_fit_melt_curve(RefSeq, cluster_df, min_fluor_subset, max_fluor_subset, condition_list, n_bootstraps=100, test=False):
    '''
    Wrapper to access clusters corresponding to `RefSeq` sequence from `cluster_df`.
    '''

    df_signal = cluster_df.loc[cluster_df.RefSeq==RefSeq][cluster_df.whichThreePrime==0]
    min_signal = min_fluor_subset.loc[min_fluor_subset.whichThreePrime==0]
    max_signal = max_fluor_subset.loc[max_fluor_subset.whichThreePrime==0]

    print('lens', len(min_signal),len(max_signal),)

    dH_list, dS_list, dG_list, Tm_list, rmse_list, values_list, signal_list = [],[],[],[],[],[],[]

    ### Bootstrap two-state fit

    for bs_ind in range(n_bootstraps):

        values = []
        signal=[]

        minFluor_df = resample_df(min_signal)
        maxFluor_df = resample_df(max_signal)
        variantFluor_df = resample_df(df_signal)

        for c in condition_list: # for each temperature

            minFluorescence = np.nanmedian(minFluor_df[c].values)
            maxFluorescence = np.nanmedian(maxFluor_df[c].values)
            variantFluorescence = np.nanmedian(variantFluor_df[c].values)

            VariantProbUnfolded = (variantFluorescence - minFluorescence)/(maxFluorescence - minFluorescence)

            values.append(VariantProbUnfolded)
            signal.append(variantFluorescence)

        if any(np.isfinite(values)):

            dH, dS, dG, Tm, rmse = fit_melt_curve(values)

        else:
            dH, dS, dG, Tm, rmse = np.nan, np.nan, np.nan, np.nan, np.nan

        dH_list.append(dH)
        dS_list.append(dS)
        dG_list.append(dG)
        Tm_list.append(Tm)
        rmse_list.append(rmse)
        values_list.append(values)
        signal_list.append(signal)

    values_list = np.array(values_list)
    signal_list = np.array(signal_list)

    dH_estimate = np.median(dH_list)
    Tm_estimate = np.median(Tm_list)

    # if test:
    #     np.savetxt("%s_example_melt_curves_%s.npy" % (RefSeq, todaysdate),values_list, )

    # non-bootstrapped cluster frac. unfolded to perform chisquare test

    T_celsius=[15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,]
    T_kelvin=[x+273.15 for x in T_celsius]

    chisquared_all_clusters = 0
    chisquared_of_median = 0

    if test: 
        print('RefSeq')
        print('Tm est: %.2f' % Tm_estimate)
        print('dH est: %.2f' % dH_estimate)

    n_observations = 0

    for T, c in list(zip(T_kelvin, condition_list)):
        c_min = np.nanmedian(min_signal[c].values)
        c_max = np.nanmedian(max_signal[c].values)

        df_signal[c+'_FracUnfolded'] = (df_signal[c] - c_min)/(c_max-c_min)

        pred_frac_unfolded = predict_two_state(dH_estimate, Tm_estimate, T)
        n_observations += (len(df_signal[c+'_FracUnfolded']) - len(df_signal[df_signal[c+'_FracUnfolded'].isna()]))
        chisquared_all_clusters += np.sum((df_signal[c+'_FracUnfolded'].values - pred_frac_unfolded)**2 / np.nanvar(df_signal[c+'_FracUnfolded'].values))
        chisquared_of_median += np.sum((np.nanmean(df_signal[c+'_FracUnfolded'].values) - pred_frac_unfolded)**2 / np.nanvar(df_signal[c+'_FracUnfolded'].values))

    chisquared_all_clusters /= (n_observations - 2)
    chisquared_of_median /= 17

    if test:
        print('Chisquared all clusters: %.2f' % chisquared_all_clusters)
        print('Chisquared of median: %.2f' % chisquared_of_median)

    # Prepare output

    output_dct={'RefSeq': RefSeq,
                'n_clusters_signal': len(df_signal),
                'dH': np.median(dH_list),
                'dH_err': np.std(dH_list),
                'dS': np.median(dS_list),
                'dS_err': np.std(dS_list),
                'dG_37C': np.median(dG_list),
                'dG_37C_err': np.std(dG_list),
                'RMSE': np.mean(rmse_list),
                'chisquared_per_dof': chisquared_of_median,
                'chisquared_all': chisquared_all_clusters,
                'Tm': np.median(Tm_list),
                'Tm_err': np.std(Tm_list),
                }

    for i, c in enumerate(condition_list):
        output_dct.update({c+'_median': np.nanmedian(values_list[:,i])})
        output_dct.update({c+'_err': np.nanstd(values_list[:,i])})

        output_dct.update({c+'_signal_median': np.nanmedian(signal_list[:,i])})
        output_dct.update({c+'_signal_err': np.nanstd(signal_list[:,i])})

    return pd.Series(output_dct)

def predict_two_state(dH, Tm, T):
    function = lambda dH, Tm, x: 1 / (1 + np.exp(dH/0.00198*(Tm**-1 - T**-1)))
    pred_fit = function(dH, Tm, T)
    return pred_fit
    

def fit_melt_curve(orig_vals, K=0.00198, epsilon = 1e-13):

    x_temps=[x+273.15 for x in [15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,]]

    #takes in a series of values, returns fit dH, dS, dG at 37 ˚C, Melting temperature (Tm)
    
    # p[0] is dH / k, p[1] is T_m
    sigmoid_f = lambda p, x, y:  (abs(y - 1 / (1 + np.exp(p[0]*(p[1]**-1 - x))))).sum()
    function = lambda p, x: 1 / (1 + np.exp(p[0]*(p[1]**-1 - x)))

    orig_vals = np.array(orig_vals)
    assert len(orig_vals) == len(x_temps)
    finite_inds = np.isfinite(orig_vals)
    finite_vals = orig_vals[finite_inds]
    vals_transform = np.clip(finite_vals, epsilon, 1-epsilon)

    x_temps = np.array(x_temps)[finite_inds]
    x_inv = [1/x for x in x_temps]

    p0 = [-100.0, 300.0] # initial values for dH, Tm

    xopt = opt.fmin(func=sigmoid_f, x0=p0, args=(x_inv, vals_transform), disp=False) # downhill simplex algorithm
    fit_dH_k, fit_Tm = xopt
    pred_fit = function(xopt, x_inv)

    #fit_Tm = 1/fit_Tm_inv
    fit_dH = fit_dH_k*K
    fit_dS = fit_dH/fit_Tm
    fit_dG_37C = fit_dH - fit_dS * (273.15+37)

    fit_RMSE = np.sqrt(np.mean(np.square(finite_vals - pred_fit)))

    return fit_dH, fit_dS, fit_dG_37C, fit_Tm, fit_RMSE

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Fit melt curve to nucleic acid hairpin library variants.

    H Wayment-Steele, E Sharma, Aug 2021

    Input: Dataframe of single clusters, with raw intensity at each temperature point (written by ProcessData.py)
    Output: .json file containing statistics bootstrapped over clusters for each variant

    """)

    p.add_argument("infile", action="store", help="Input json or path to jsons.") # each row has a RefSeq, and one field per temperature 
    p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
    p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
    p.add_argument('--min_fluor_seqs', action='store', help='List of constructs to use to set minimum fluorescence.')
    p.add_argument('--max_fluor_seqs', action='store', help='List of constructs to use to set max fluorescence.')

    p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_out.json')
    p.add_argument('--n_bootstraps', action="store", type=int, default=100, help='Number of bootstraps for resampling each variant.')

    args = p.parse_args()
    print(args)

    if args.parallel:
        pandarallel.initialize()

    if args.outfile:
        out_prefix = args.outfile
    else:
        out_prefix = os.path.basename(args.infile).replace('.json','')

#### Read in ####

    if args.min_fluor_seqs is None:
        min_fluor_seqs = ['GGCCGCGGAAACGCGGCC']
    else:
        min_fluor_seqs=[]
        with open(args.min_fluor_seqs,'r') as f:
            for lin in f.readlines():
                min_fluor_seqs.append(lin.strip())

    if args.max_fluor_seqs is None:
        max_fluor_seqs = ['A'*18]
    else:
        max_fluor_seqs=[]
        with open(args.max_fluor_seqs,'r') as f:
            for lin in f.readlines():
                max_fluor_seqs.append(lin.strip())

    if args.infile.endswith('.json') or args.infile.endswith('.json.zip'):
        df = pd.read_json(args.infile)
    elif args.infile.endswith('.pkl'):
        df = pd.read_pickle(args.infile)
    else:
        #assume path to directory
        files = glob(args.infile+'/*.json.zip')

        if args.test:
            files = files[:1]

        print('found files to read in:')
        for x in files:
            print(x)

        df = pd.DataFrame()
        for x in files:
            tmp = pd.read_json(x)
            df = df.append(tmp, ignore_index=True)

    print('Read in %s' % args.infile)

### Hardcoded! ####
    green_conditions = [x for x in list(df.keys())[:-2] if 'Green' in x]
    red_conditions = [x for x in list(df.keys())[:-2] if 'Red' in x]

    #normalize each cluster to its Red channel
    green_conditions_to_normalize = green_conditions[4:-1]
    green_conditions_norm = [x+'_redNorm' for x in green_conditions_to_normalize]

#### Normalize ####

    df = df.dropna(subset=red_conditions)
    print("Removed clusters with NaNs in red")

    for cond in red_conditions:
        df[cond] = clip_normalize(df[cond])

    print('Removed outliers from red')

    for cond in green_conditions_to_normalize:
        corresponding_red_condition = cond.replace('Green','Red')
        df[cond+'_redNorm'] = df[cond] / df[corresponding_red_condition]

    norm_conditions = [x+'_redNorm' for x in green_conditions_to_normalize]

    max_fluor_subset = df.loc[df.RefSeq.isin(max_fluor_seqs)]
    min_fluor_subset = df.loc[df.RefSeq.isin(min_fluor_seqs)]

    for cond in norm_conditions:
        min_fluor_subset[cond] = clip_normalize(min_fluor_subset[cond])
        max_fluor_subset[cond] = clip_normalize(max_fluor_subset[cond])

    min_fluor_subset.to_json('min_fluor_subset.json.zip')
    max_fluor_subset.to_json('max_fluor_subset.json.zip')

    ### Get list of unique variants ###

    unique_RefSeqs = list(df['RefSeq'].unique())

    if args.test: unique_RefSeqs = unique_RefSeqs[:20]

    RefSeq_list = pd.DataFrame({'RefSeq': unique_RefSeqs})


    if args.parallel:
        output_df = RefSeq_list.parallel_apply(lambda row: bootstrap_fit_melt_curve(row['RefSeq'], df, min_fluor_subset,max_fluor_subset, norm_conditions, n_bootstraps = args.n_bootstraps, test=args.test), axis=1, result_type='expand')
    else:
        output_df = RefSeq_list.progress_apply(lambda row: bootstrap_fit_melt_curve(row['RefSeq'], df, min_fluor_subset,max_fluor_subset, norm_conditions, n_bootstraps = args.n_bootstraps, test=args.test), axis=1, result_type='expand')

    outfile = out_prefix+ '_BootstrappedFits_' + todaysdate + '.json.zip'
    output_df.to_json(outfile)
    print('Saved Bootstrapped Fits to %s' % outfile)

