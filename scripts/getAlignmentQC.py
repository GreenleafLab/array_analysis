"""
Functions and scripts for QC plots and tables of the MiSeq sequence alignment result.
Yuxi Ke, Nov 2021
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
sns.set_context('paper')


class FileNameMaker:
    def __init__(self, fig_dir, experiment_name):
        self.fig_dir = fig_dir
        self.experiment_name = experiment_name

    def get_plot_path(self, short_name, ext='.png'):
        return os.path.join(self.fig_dir, self.experiment_name + '_' + short_name + ext)


def get_stats(CPseq_file, stats_file):
    df = pd.read_csv(CPseq_file, delimiter='\t')
    df['GC content'] = df.apply(lambda row: (row['RefSeq'].count('G')+row['RefSeq'].count('C'))/len(row['RefSeq']), axis=1) 
    tmp = df.groupby(['RefSeq', 'Series']).size().rename(columns={'0':'n_cluster'})

    tmp.to_csv(stats_file)


def get_num_per_series(df, annotation):
    """
    TODO: Realign with the new version of annotation file to resolve name differences
    """
    number_per_series = pd.DataFrame(df.groupby('series')['RefSeq'].count()).merge(annotation.groupby('series').count(), on='series').iloc[:,:2].rename(columns={'clusterID': 'cluster_number', 'RefSeq': 'variant_number'})
    number_per_series['cluster-variant ratio'] = number_per_series['cluster_number'] / number_per_series['variant_number']


def is_present_over_cutoff(refseq, alignstats, cutoff=20):
    if not (refseq in alignstats['RefSeq'].tolist()):
        #print('not there')
        return False
    else:
        if alignstats[alignstats['RefSeq'] == refseq].iloc[0,:]['n_cluster'] <= cutoff:
            return False
        else:
            return True


def plot_num_clusters_per_variant(align_stats, plot_path):
    fig, ax = plt.subplots(figsize=(6,14))
    plt.title('#clusters per variant')
    sns.boxplot(data=alignstats, y='Series', x='n_cluster', orient='h', palette='Pastel1', width=0.4)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)


def plot_percent_variant_presented(coverage, plot_path):
    """
    Args:
        coverage - df, columns = ['percentage', 'percentage_pass_filter']
    """
    fig, ax = plt.subplots(1, 2, figsize=(6,9), sharey=True)
    sns.barplot(y=coverage.index, x=coverage.percentage, palette='tab10', ax=ax[0])
    plt.xlabel('%')
    ax[0].set_title('% variant presented')
    sns.barplot(y=coverage.index, x=coverage.percentage_pass_filter, palette='tab10', ax=ax[1])
    ax[1].set_title('% variant presented > 20 clusters')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)


def process_stats_file(stats_file, fn_maker):
    print('Reading align stats...\n')
    alignstats = pd.read_csv(stats_file)
    if '0' in alignstats.columns:
        alignstats = alignstats.rename(columns={'0':'n_cluster'})
    #alignstats = pd.DataFrame(alignstats.groupby('RefSeq').sum()['n_cluster']).merge(annotation, how='left', on='RefSeq')
    print('Number of unique variants on the chip: %d\nNumber of clusters: %d' % (len(np.unique(alignstats['RefSeq'])), alignstats['n_cluster'].sum()))

    print('\n\nPlotting number of clusters per variant')
    plot_num_clusters_per_variant(alignstats, plot_path=fn_maker.get_plot_path('num_cluster_per_variant', ext='.pdf'))
    plot_num_clusters_per_variant(alignstats, plot_path=fn_maker.get_plot_path('num_cluster_per_variant', ext='.png'))


def process_CPseq_file(CPseq_file, annotation_file, fn_maker)

if __name__ == "__main__":
    if 'snakemake' in globals():
        experiment_name = snakemake.config['experimentName']
        stats_file = snakemake.input['align_stats']
        CPseq_file = snakemake.input['aligned_CPseq']
        annotation_file = snakemake.input['annotation_file']
        fig_dir = snakemake.output['fig_dir']
    else:
        # test
        datadir = r'/scratch/groups/wjg/kyx/NNNlib2b_Nov11'
        experiment_name = 'NNNlib2b'
        stats_file = '/scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/aligned/ConsensusReads_20211118_exact_STATS.csv'
        CPseq_file = '/scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/aligned/ConsensusReads_20211118_exact.CPseq'
        annotation_file = '/scratch/groups/wjg/kyx/NNNlib2b_Nov11/data/reference/NNNlib2b_annotation_nupack.tsv'
        fig_dir = '/scratch/groups/wjg/kyx/NNNlib2b_Nov11/fig/align_QC'

    fn_maker = FileNameMaker(fig_dir, experiment_name)

    process_stats_file(stats_file, fn_maker)

    #df = pd.read_csv(CPseq_file, sep='\t')
