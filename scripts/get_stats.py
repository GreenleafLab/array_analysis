import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], delimiter='\t')
df['GC content'] = df.apply(lambda row: (row['RefSeq'].count('G')+row['RefSeq'].count('C'))/len(row['RefSeq']), axis=1) 
tmp = df.groupby(['RefSeq','whichThreePrime','series']).size()

tmp.to_csv('%s_STATS.csv' % sys.argv[1].replace('.csv',''))
