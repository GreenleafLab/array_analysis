import nwalign3 as nw
import os, subprocess, argparse
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from tqdm import tqdm

tqdm.pandas()

pandarallel.initialize(progress_bar=True)
# gapPenalty1 = 0
# gapExtension1 = 0
# gapPenalty2 = -1
# gapExtension2 = 0
# scoringMatrix = "data/reference/NUC.4.4"

def getScorePvalue(nwScore, m, n, k=0.0097, l=0.5735, nwScoreScale=0.2773):
    """Get pvalue of extreme value distribution which model's likelihood of achieveng a more extreme
    alignment score for two sequences of length m and n. K and l are empirically determined.
    nwScoreScale is a factor applied to scores of NUC4.4 matrices (from MATLAB nwalign)"""
    u = np.log(k*m*n)/l
    score = nwScore*nwScoreScale
    return (1 - np.exp(-np.exp(-l*(score-u))))

def calcMinQScore(libRegion):
    q_scores = [ord(x)-33 for x in libRegion]
    return np.min(q_scores)

def getAlignment(seq1_in, seq2_in, gapPenalty, gapExtension):
    seq1_out, seq2_out = nw.global_align(seq1_in, seq2_in, gap_open=gapPenalty, gap_extend=gapExtension, matrix=scoringMatrix)
    score = nw.score_alignment(seq1_out, seq2_out, gap_open=gapPenalty,  gap_extend=gapExtension, matrix=scoringMatrix)
    pvalue = getScorePvalue(score, len(seq1_in), len(seq2_in))
    return seq1_out, seq2_out, score, pvalue

def get_hairpin(row, fiveprime_region ='GCTGTTGAAGGCTCGCGATGCACACGCTCTGGTACAAGGAA',
                        threeprime_region='AAGGCACTGGGCAATACGAGCTCAAGCCAGTCTCGCAGTCC', 
                        debug=False, pValueCutoff=1e-3):

    '''
    Align read to two flanking regions and return region in between them.

    Input:
        row - row of dataframe containing fields for `sequence` and `phred`
        fiveprime_region - str, 5' region containing Cy3 binding site and 2 linker 'A'
        threeprime_region - str, either 2 alternatives (historical design, one with quench and one without)
            or only 1 (current disign, quench binding site and extra 'A's). Put it in a list to keep consistant
    Output:
        libregion: detected region between two flanking regions (str)
        phred_libregion: phred of detected region

    If sequence fails to align: returns NaN's for all
    '''

    sequence=row['sequence']
    seq1a, seq1b, score1, pvalue1 = getAlignment(sequence, fiveprime_region, 0, 0)

    if pvalue1 < pValueCutoff:

        seq2a, seq2b, score2, pvalue2 = getAlignment(sequence, threeprime_region, -1, 0)

        if pvalue2 < pValueCutoff:
            lib_start = seq1b.find('-')
            lib_end = seq2b.find(threeprime_region[0])

            libregion = sequence[len(fiveprime_region):lib_end]
            phred_libregion = row['phred'][len(fiveprime_region):lib_end]

            return libregion, phred_libregion
        else:
            return np.nan, np.nan

    else:
        return np.nan, np.nan
    

def getLibraryRef(ex_seq, sort_lib_seqs, beam=100, second_try = False, exact=False,
                  gapPenalty1 = -1, gapExtension1 = -1, pValueCutoff = 1e-6, debug=False):
    """
    Returns:
        RefSeq
        pValue
    """
    if ex_seq=='':
        return ''
    
    idx = np.searchsorted(sort_lib_seqs, ex_seq)

    if idx < len(sort_lib_seqs) and sort_lib_seqs[idx] == ex_seq:
        return sort_lib_seqs[idx], 0

    else:
        if exact:
            return np.nan, np.nan # only returning if exact match
        else:

            min_index = max(0, idx-beam)
            max_index = min(len(sort_lib_seqs), idx+beam)

            scores, pvalues = [],[]

            for trial_seq in sort_lib_seqs[min_index:max_index]:
                seq1a, seq1b = nw.global_align(ex_seq, trial_seq, gap_open=gapPenalty1, gap_extend=gapExtension1, matrix=scoringMatrix)
                score1 = nw.score_alignment(seq1a, seq1b, gap_open=gapPenalty1, gap_extend=gapExtension1, matrix=scoringMatrix)
                pvalue1 = getScorePvalue(score1, len(ex_seq), len(trial_seq))
                if debug:
                    print(seq1a)
                    print(seq1b)
                    print(score1, pvalue1)

                scores.append(score1)
                pvalues.append(pvalue1)

            winner = np.argmax(scores)

            if pvalues[winner] < pValueCutoff:
                return sort_lib_seqs[min_index:max_index][winner], pvalues[winner]

            #else:
            #    if not second_try:
            #        return getLibraryRef(ex_seq, sort_lib_seqs, beam=100000, second_try = True)
            else:
                return np.nan, np.nan # didn't find a matching sequence

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("cpseq", help="CPseq file containing paired-end reads.")
    parser.add_argument("--library", action='store',  help="CSV file containing library information. Must have column called 'RefSeq'.")
    parser.add_argument("--exact", action='store_true',  help="Only return sequences that are exact match to library RefSeqs.")
    parser.add_argument("--beam", action='store', default=1000, help='beam search width.')
    parser.add_argument("--OligoPValue", action='store', default=1e-3, help='P-value cutoff for aligning reads to fluor and quench oligo.')
    parser.add_argument("--LibPValue", action='store', default=1e-6, help='P-value cutoff for aligning libregions to library.')
    parser.add_argument('-o', action='store',help='name of output CPseq')
    parser.add_argument('--scoringMatrix', action='store', default='NUC.4.4', help='path to the scoring matrix file')
    parser.add_argument('--fiveprimeRegion', action='store', default='GCTGTTGAAGGCTCGCGATGCACACGCTCTGGTACAAGGAA')
    parser.add_argument('--threeprimeRegion', action='store', default='AAGGCACTGGGCAATACGAGCTCAAGCCAGTCTCGCAGTCC')
    

    args = parser.parse_args()

    if args.o is None:
        args.o = 'annotated_output.csv'
    global scoringMatrix
    scoringMatrix = args.scoringMatrix
    print('Using scoring matrix %s' % scoringMatrix)
    exact = args.exact
    if args.beam <= 1:
        exact = True

    library = pd.read_csv(args.library)
    #clean RefSeq column
    library['RefSeq'] = [x.upper().replace('U','T') for x in library['RefSeq']]
    # ensure ref seqs are sorted
    sort_lib_seqs = list(np.sort(list(library['RefSeq'])))

    print('Read library')

    if os.path.exists('%s_libregion_only.json.zip' % args.o):
        df = pd.read_json('%s_libregion_only.json.zip' % args.o)

    else:
        df = pd.read_csv(args.cpseq, names=['clusterID','sequence','phred'], delimiter='\t')
        df['clusterID'] = [x.replace('@','').replace(' 1:N:0:1','') for x in df['clusterID']]
        print('Read CPseq in')
        print(df.head())
        
        print("Extracting variable region from between Cy3' and quench' regions....")
        df[['libRegion', 'libRegionPhred']] = df.apply(lambda row: get_hairpin(row, pValueCutoff=args.OligoPValue, fiveprime_region=args.fiveprimeRegion, threeprime_region=args.threeprimeRegion), axis=1, result_type='expand')

        df.to_json('%s_libregion_only.json.zip' % args.o)
        #df.iloc[:1000].to_csv('%s_libregion_only_test.csv' % args.o,index=False)

    n_passed_alignment = len(df.loc[~df['libRegion'].isna()])

    print("%d/%d (%.2f) passed aligning to Cy3' quench' regions" % (n_passed_alignment, len(df), 100*n_passed_alignment/len(df)))
    
    df = df.loc[~df['libRegion'].isna()]
    df['len_libRegion'] = df.apply(lambda row: len(row['libRegion']), axis=1)
    df = df.loc[df['len_libRegion']>0]
    df['minQScore'] = df.parallel_apply(lambda row: calcMinQScore(row['libRegionPhred']), axis=1)

    print('Len prior to minQscore filtering: ', len(df))
    df = df.loc[df['minQScore']>=20]
    print('Length with minQscore>20', len(df))

    uniqueLibRegions = pd.DataFrame()
    uniqueLibRegions['libRegion'] = df['libRegion'].unique()

    print("Found %d unique libRegions to align" % len(uniqueLibRegions))

    # match each libRegion to the most likely reference sequence from the library
    print("Aligning library region to reference sequences and matching to most likely ref sequence ....")
    uniqueLibRegions[['RefSeq','RefSeqPValue']] = uniqueLibRegions.parallel_apply(lambda row: getLibraryRef(row['libRegion'],sort_lib_seqs, exact=exact, pValueCutoff=args.LibPValue, beam=args.beam), axis=1, result_type='expand')

    #merge the library data to the unique_libregion data
    uniqueLibRegions = uniqueLibRegions.merge(library, on='RefSeq')

    #merge the unique_libregion data to the original CPseq data
    df = df.merge(uniqueLibRegions, on='libRegion')

    n_passed_alignment2 = len(df.loc[~df['RefSeq'].isna()])

    print("%d/%d (%.2f) passed aligning to library RefSeqs" % (n_passed_alignment2, n_passed_alignment, n_passed_alignment2/n_passed_alignment))

    # remove libRegions that did not align to a ref seq
    df= df.loc[~df['RefSeq'].isna()]

    df.to_csv(args.o, sep='\t', index=False)