"""
Functions and scripts for QC plots and tables of the alignment result.
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

    def get_plot_path(self, short_name, ext='.pdf'):
        return os.path.join(self.fig_dir, self.experiment_name + '_' + short_name + ext)


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
    fig, ax = plt.subplots(figsize=(14,4))
    plt.title('#clusters per variant')
    sns.boxplot(data=alignstats, x='series', y='n_cluster', palette='Pastel1', width=0.4)
    plt.ylim([0, 1000])
    plt.xticks(rotation=45)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    if 'snakemake' in globals():
        experiment_name = snakemake.config['experimentName']
        stats_file = snakemake.input['align_stats']
        CPseq_file = snakemake.input['aligned_CPseq']
        refseq_file = snakemake.input['refseq_file']
        fig_dir = snakemake.output['fig_dir']
    else:
        # test
        experiment_name = 'NNNlib2b'
        stats_file = '~/workspace/nnn/data/NNNlib2b/ConsensusReads_20211115_exact_STATS.csv'
        CPseq_file = '~/workspace/nnn/data/NNNlib2b/ConsensusReads_20211115_exact.CPseq'
        refseq_file = '~/workspace/nnn/data/NNNlib2b/NNNlib2b_RefSeqs.csv'
        annotation_file = '~/workspace/nnn/data/NNNlib2b_annotation_nupack.txt'
        fig_dir = '/Users/yuxi/workspace/nnn/fig/align_QC'

    fn_maker = FileNameMaker(fig_dir, experiment_name)

    print('\n\nReading library annotations')
    annotation = pd.read_csv(annotation_file, sep='\t')
    annotation['RefSeq'] = annotation['RefSeq'].apply(lambda s: s.upper().replace('U', 'T'))

    print('Reading align stats...')
    alignstats = pd.read_csv(stats_file)
    alignstats = alignstats.rename(columns={'0':'n_cluster'})
    alignstats = pd.DataFrame(alignstats.groupby('RefSeq').sum()['n_cluster']).merge(annotation, how='left', on='RefSeq')
    print('Number of unique variants on the chip: %d\nNumber of clusters: %d' % (len(np.unique(alignstats['RefSeq'])), alignstats['n_cluster'].sum()))


    print('\n\nPlotting number of clusters per variant')
    plot_num_clusters_per_variant(alignstats, plot_path=fn_maker.get_plot_path('num_cluster_per_variant'))


    df = pd.read_csv(CPseq_file, sep='\t')
