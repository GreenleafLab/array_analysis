import os
import pandas as pd

def fill_sample_sheet(mapfile):
    """
    Fill tabular csv sample sheet.
    Args:
        mapfile - str
            a csv file with the column 'image_set' 
            that's directory names e.g. '15_Green_15C'
    """
    df = pd.read_csv(mapfile)
    df['channel'] = df['image_set'].apply(lambda s: s.split('_')[1])
    df['condition'] = df['image_set'].apply(lambda s: s.split('_')[2])
    df['is_temperature'] = df['condition'].apply(lambda s: s.endswith('C'))

    df.to_csv(mapfile)



def parse_mapfile(mapfile):
    """
    Args:
        mapfile - str
    Returns:
        fluordir - str
        seriesdir - str
        dirname - List[str], all fluor dirs
        condition - List[str]
    """
    with open(mapfile, 'r') as mf:
        mapinfo = mf.readlines()

    mapinfo = [l.strip('\n') for l in mapinfo]

    fluordir = mapinfo[0]
    seriesdir = mapinfo[1]
    # channel = [s.split("_")[1] for s in mapinfo[1:]]
    dirname = mapinfo[2:]
    condition = [s.split("_")[-1].strip('\n') for s in dirname]

    return fluordir, seriesdir, condition, dirname

def parse_fluorfiles_from_mapfile(mapfile):
    """
    Returns:
        fluorfile - List[str]
        seriesfile - List[str]
    """
    
    fluordir, seriesdir, condition, dirname = parse_mapfile(mapfile)
    fluorfolder = [os.path.join(fluordir, s) for s in dirname]
    
    fluorfile, seriesfile = [], []
    for i,dir in enumerate(fluorfolder):
        fluorfile.extend([os.path.join(dir, s) for s in os.listdir(dir) if s.endswith('.CPfluor')])
        seriesfile.extend([os.path.join(seriesdir, dirname[i], s.strip('.CPfluor')+'.CPseries') for s in os.listdir(dir) if s.endswith('.CPfluor')])

    return fluorfile, seriesfile

fill_sample_sheet(r'config/nnnlib2_map.csv')
#print(parse_mapfile(r'../config/test.map'))
#print(parse_fluorfiles_from_mapfile(r'../config/test.map'))
