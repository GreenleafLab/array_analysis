"""
Normalize pooled CPseries files for NNNlib2b.
Input:
    Merged CPseries file
    A mapfile csv file that instructs which columns to normalize and keep
    An annotation tsv file that determines which variants to use as initial Fmax and Fmin controls
Output:
    Normalized data with index 'clusterID' and columns final_norm_conditions
    where sequence number and temperatures are in the `final_norm_conditions` header

Yuxi Ke, Nov 2021
"""
from typing_extensions import final
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
sns.set_context('paper')

import os
import warnings


def clip_normalize(data):
    '''
    Used when normalizing signal.
    '''

    min_fluor = np.nanpercentile(data,1)
    max_fluor = np.nanpercentile(data,99)
    return np.clip(data, min_fluor, max_fluor)


def get_conditions_from_mapfile(mapfile, channel):
    """
    Args:
        mapfile - str, file name
        channel - str, {'green', 'red'}
    Return:
        condition_list - List[str]
    """
    metadata = pd.read_csv(mapfile)
    return metadata[metadata['curve_channel'] == channel]['condition'].tolist()


def get_xdata_from_condition(condition):
    """
    Args:
        condition - str, e.g. Green13_25
    Returns:
        xdata - str, temperature in Kalvin
    """
    xdata = float(condition.split('_')[1]) + 273.15
    return '%.2f'%xdata


def get_long_and_stem_refseq(annotation):

    control_refseq = np.unique(annotation[annotation['ConstructType'] == 'Control']['RefSeq'])
    long_refseq = [seq for seq in control_refseq if len(seq) >= 39]
    stem_refseq = np.unique(annotation[annotation['ConstructType'] == 'Super_Stable_Stem']['RefSeq'])

    return long_refseq, stem_refseq


def get_refseq_median(refseq, df, conditions):

    warnings.filterwarnings("ignore")

    temperature = [cond.split('_')[1] for cond in conditions]
    vardf = df[df['RefSeq'] == refseq]
    refseq_median = np.nanmedian(vardf[conditions].values, axis=0)
    
    return temperature, refseq_median


def get_control_refseq_medians_and_plot(clean_df, long_refseq, stem_refseq, conditions, fig_path, ylim=None):
    n_long, n_stem = len(long_refseq), len(stem_refseq)
    #print('%d long, %d stem controls' % (n_long, n_stem))

    long_median = np.zeros((n_long, len(conditions)))
    for i in range(n_long):
        _, long_median[i,:] = get_refseq_median(long_refseq[i], clean_df, conditions)

    stem_median = np.zeros((n_stem, len(conditions)))
    for i in range(n_stem):
        temperature, stem_median[i,:] = get_refseq_median(stem_refseq[i], clean_df, conditions)

    fig, ax = plt.subplots()
    plt.plot(temperature, long_median.T, 'orange')
    plt.plot(temperature, stem_median.T, 'purple')
    plt.ylim(ylim)
    plt.xlabel('Temperature (Celcius)')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    
    return long_median, stem_median


def get_good_control_refseq(control_median, control_refseq, percentile_cutoff=5, n_outlier_temp_point_thresh=1):
    """
    Args:
        control_median - (n_refseq, n_temperature)
        control_refseq - (n_refseq,) list of str
        percentile - float, double side percentile cutoff
        n_outlier_temp_point_thresh - int, how many temperature points in one curve is allowed to lie outside percentile cutoff
    Returns:
        good_control_refseq - List[str]
    """
    warnings.filterwarnings("ignore")

    good_control = (control_median > np.nanpercentile(control_median, percentile_cutoff, axis=0)).sum(axis=1) >= (control_median.shape[1] - n_outlier_temp_point_thresh)
    good_control = np.logical_and(good_control, (control_median < np.nanpercentile(control_median, 100 - percentile_cutoff, axis=0)).sum(axis=1) >= (control_median.shape[1] - n_outlier_temp_point_thresh))
    good_control_refseq = np.array(control_refseq)[good_control]

    return list(good_control_refseq)


def plot_example_refseqs(refseqs, df, conditions, annotation, cmap, ylim=None):
    fig, ax = plt.subplots()
    for refseq in refseqs:
        dG_37C = annotation[annotation['RefSeq'] == refseq]['dG_37C'].values[0]
        dG_37C = np.clip((dG_37C - (-10)) / (0 - (-10)), 0, 1)
        temperature, example_median = get_refseq_median(refseq, df, conditions)
        
        ax.plot(temperature, example_median, c=cmap(dG_37C))
    ax.set_ylim(ylim)


def plot_color_coded_WC_examples(clean_df, annotation, final_norm_conditions, fig_path, n_variant_to_plot=30):
    np.random.seed(42)
    cmap = sns.cubehelix_palette(start=1, rot=-.8, reverse=True, as_cmap=True)
    WC_refseqs = annotation[annotation['Series'] == 'WatsonCrick']['RefSeq'].values[np.random.choice((annotation['Series'] == 'WatsonCrick').sum(), n_variant_to_plot)]
    plot_example_refseqs(WC_refseqs, clean_df, final_norm_conditions, annotation, cmap=cmap, ylim=[-.2, 1.2])
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')


if __name__ == '__main__':

    if 'snakemake' in globals():
        # Running from snakemake
        CPseries_file = snakemake.input['CPseries_file']
        mapfile = snakemake.input['mapfile']
        annotation_file = snakemake.input['annotation']
        green_norm_condition = snakemake.params['green_norm_condition']
        figdir = snakemake.params['figdir']
        out_file = snakemake.output['out_file']
        xdata_file = snakemake.output['xdata_file']
        ext = snakemake.params['ext']
        variant_col = snakemake.params['variant_col']
    else:
        # Runnig as a test
        CPseries_file = r'C:\Users\Yuxi\workspace\NNNlib2b_CPseries\NNNlib2b_DNA_20211022.pkl'
        mapfile = r'C:\Users\Yuxi\workspace\NNNlib2b_CPseries\nnnlib2b_map.csv'
        annotation_file = r'C:\Users\Yuxi\workspace\NNNlib2b_CPseries\NNNlib2b_annotation_nupack.txt'
        green_norm_condition = 'Green07_PostCy3'
        figdir = r'C:\Users\Yuxi\workspace\NNNlib2b_CPseries\fig\20211123'
        out_file = r'C:\Users\Yuxi\workspace\NNNlib2b_CPseries\NNNlib2b_DNA_20211022_normalized.pkl'
        xdata_file = r'C:\Users\Yuxi\workspace\NNNlib2b_CPseries\NNNlib2b_DNA_20211022_xdata.txt'
        ext = '.png'
        variant_col = 'SEQID'

    # Load the data and condition names
    clean_df = pd.read_pickle(CPseries_file).dropna()#.dropna(axis=0, thresh=5)
    metadata = pd.read_csv(mapfile)
    annotation = pd.read_csv(annotation_file, sep='\t')
    green_conditions = get_conditions_from_mapfile(mapfile, 'green')
    print('Data loaded.')

    # create figdir
    if not os.path.isdir(figdir):
        os.makedirs(figdir)

    # long control and super stem refseqs
    long_refseq, stem_refseq = get_long_and_stem_refseq(annotation)

    # GreenNorm
    green_norm_conditions = [x+'_greenNorm' for x in green_conditions]
    for cond in green_conditions:
        clean_df[cond + '_greenNorm'] = clean_df[cond] / clip_normalize(clean_df[green_norm_condition])

    # Calculate medians and plot after green normalization
    fig_path = os.path.join(figdir, 'green_normalized_curves_01-all_controls' + ext)
    norm_long_median, norm_stem_median = get_control_refseq_medians_and_plot(clean_df, long_refseq, stem_refseq, green_norm_conditions, fig_path=fig_path, ylim=[0, 3])
    
    # Pick the good max and min control variants
    good_long_refseq = get_good_control_refseq(norm_long_median, long_refseq, percentile_cutoff=15, n_outlier_temp_point_thresh=1)
    good_stem_refseq = get_good_control_refseq(norm_stem_median, stem_refseq, percentile_cutoff=1, n_outlier_temp_point_thresh=15)
    print('%d/%d selected for long Fmax controls' % (len(good_long_refseq), len(long_refseq)))
    print('%d/%d selected for stem Fmin controls' % (len(good_stem_refseq), len(stem_refseq)))

    # Calculate initial Fmax and Fmin
    fig_path = os.path.join(figdir, 'green_normalized_curves_02-only_selected_good_controls' + ext)
    max_median, min_median = get_control_refseq_medians_and_plot(clean_df, good_long_refseq, good_stem_refseq, green_norm_conditions, fig_path=fig_path, ylim=[0, 3])
    max_median = np.median(max_median, axis=0)
    min_median = np.median(min_median, axis=0)
    max_min_median = max_median - min_median
    print('Fmax and Fmin calculated.')

    # Normalize green normed signal to calculated initial Fmax and Fmin
    final_norm_conditions = [x+'_norm' for x in green_conditions]
    for i,cond in enumerate(green_conditions):
        clean_df[cond + '_norm'] = (clean_df[cond + '_greenNorm'] - min_median[i]) / max_min_median[i]

    # Sanity check plots
    fig_path = os.path.join(figdir, 'max_min_normalized_curves_01-good_controls' + ext)
    _,_ = get_control_refseq_medians_and_plot(clean_df, long_refseq=good_long_refseq, stem_refseq=good_stem_refseq, conditions=final_norm_conditions, fig_path=fig_path, ylim=[-.2, 1.5])
    fig_path = os.path.join(figdir, 'max_min_normalized_curves_02-example_WC_melt_curves' + ext)
    plot_color_coded_WC_examples(clean_df, annotation, final_norm_conditions, fig_path, n_variant_to_plot=40)

    # Save normalized data
    clean_df[['clusterID', variant_col] + final_norm_conditions].set_index('clusterID').to_pickle(out_file)
    # clean_df[['clusterID'] + final_norm_conditions].set_index('clusterID').to_pickle(out_file)
    xdata = [get_xdata_from_condition(s) + '\n' for s in green_conditions]
    with open(xdata_file, 'w+') as fh:
        fh.writelines(xdata)
    print('Saved normalized data to %s' % out_file)
    print('Saved xdata to %s' % xdata_file)
