import numpy as np
import pandas as pd

import seaborn as sns
sns.set_style('ticks')
sns.set_context('paper')

import matplotlib.pyplot as plt
from glob import glob

with open('conditions_list','r') as f:
    conditions = [x.strip() for x in f.readlines()]

global_data = pd.DataFrame()

for condition in conditions:
    green_imgs = glob('%s/*.tif' % condition)
    red_imgs = glob('%s/*.tif' % condition.replace('Green','Red'))

    for img in green_imgs:
        I = plt.imread(img).flatten()
        tile_id = int(img.split('_tile')[1].split('_')[0])
        global_data = global_data.append({'Condition':condition, 'TileID': tile_id,'Color':'Green', 'MedianIntensity':np.median(I),'StdDevIntensity': np.std(I)}, ignore_index=True)

    for img in red_imgs:
        I = plt.imread(img).flatten()
        tile_id = int(img.split('_tile')[1].split('_')[0])
        global_data = global_data.append({'Condition':condition.replace('Green','Red'), 'TileID': tile_id,'Color':'Red', 'MedianIntensity':np.median(I),'StdDevIntensity': np.std(I)}, ignore_index=True)

    print("Processed", condition)


global_data.to_csv('quicklook_Tiles.csv',index=False)

plt.figure(figsize=(12,8))
plt.subplot(1,2,1)
sns.stripplot(y='Condition', x='MedianIntensity', data=global_data.loc[global_data.Color=='Green'])
plt.title('Green Median Intensity')

plt.subplot(1,2,2)
sns.stripplot(y='Condition', x='MedianIntensity', data=global_data.loc[global_data.Color=='Red'])
plt.title('Red Median Intensity')

plt.tight_layout()
plt.savefig('quicklook_Tiles.pdf',bbox_inches='tight')

