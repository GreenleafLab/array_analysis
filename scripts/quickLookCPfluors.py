import numpy as np
import pandas as pd

import seaborn as sns
sns.set_style('ticks')
sns.set_context('paper')

import matplotlib.pyplot as plt
from glob import glob

global_data = pd.DataFrame()

files = sorted(glob('*/*_tile1_*.CPfluor'))

for fil in files:
    condition = fil.split('/')[0]
    tile_id = fil.split('_tile')[1].split('_')[0]
    dat = np.loadtxt(fil, usecols=8, delimiter=':')
    if 'Green' in condition:
        color='Green'
    elif 'Red' in condition:
        color='Red'
    else:
        print('Found no color in condition')
        color='None'
    global_data = global_data.append({'Condition':condition, 'TileID': tile_id,'Color':color, 'MedianIntensity':np.median(dat),'StdDevIntensity': np.std(dat)}, ignore_index=True)
    print('Processed %s, Med Intensity %.2f' % (condition, np.median(dat)))

global_data.to_csv('quicklook_CPfluors.csv',index=False)

plt.figure(figsize=(12,8))
plt.subplot(1,2,1)
sns.stripplot(y='Condition', x='MedianIntensity', data=global_data.loc[global_data.Color=='Green'])
plt.title('Green Median Intensity')

plt.subplot(1,2,2)
sns.stripplot(y='Condition', x='MedianIntensity', data=global_data.loc[global_data.Color=='Red'])
plt.title('Red Median Intensity')

plt.tight_layout()
plt.savefig('quicklook_CPfluors.pdf',bbox_inches='tight')

