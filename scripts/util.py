import os
import pandas as pd

def get_series_tile_filenames(seriesdir):
    """
    Returns:
        tile_names - List[str]
    """
    tile_names = []
    return tile_names


def write_to_hdf5(experiment_nm, series_tile=[], h5nm=None, clean=True):
    """
    Args:
        experiment_nm - str, e.g. 'NNNlib2b_DNA'
        series_tile - List[str], list of CPseries file names
        h5nm - str, output file name
        clean - bool, whether to drop all nan and duplicate clusters
    """
    if h5nm is None:
        if not clean:
            h5nm = experiment_nm + '.h5'
        else:
            h5nm = experiment_nm + '_clean.h5'
        
    with pd.HDFStore(h5nm) as h5:
        for i in range(1,19):
            print('Reading tile %03d' % i)
            tile = pd.read_csv(series_tile[i])
            
            if clean:
                # drop clusters that only have NaNs
                # drop non-unique clusters
                tile.dropna(how='all', subset=conditions, inplace=True)
                tile.drop_duplicates(subset='clusterID', keep=False, ignore_index=True, inplace=True)
            
            if i == 1:
                h5.put('signal', tile, format='table', data_columns=True)
            else:
                h5.append('signal', tile, format='table', data_columns=True)
                
            print('Wrote tile %03d\n' % i)
                
    print('Wrote to hdf5 file %s' % h5nm)


def convert_hdf5_to_pickle(hdf5_file, pickle_file, drop_duplicate=True):

    all_df = pd.read_hdf(hdf5_file, 'signal', mode='r+')
    if drop_duplicate:
        all_df.drop_duplicates(subset='clusterID', keep=False, ignore_index=True, inplace=True)
    all_df.to_pickle(pickle_file)


def get_fluor_names_from_mapfile(mapfile, tifdir, fluordir):
    """
    Get just 1 example tile fluor file name per condition for snakemake rule `quantify_images`
    Args:
        mapfile - str
            a csv file with the column 'condition'
        tifdir - str, absolute directory for tif files, e.g. /.../data/images/
        fluordir - str, absolute directory for CPfluor dirs, e.g. /.../data/fluor/
    Returns:
        fluor_names - List[str]
    """
    df = pd.read_csv(mapfile)
    
    fluor_names = []
    for condition in df['condition']:
        tif_condition_dir = os.path.join(tifdir, condition)
        try:
            tif_fn = next(s for s in os.listdir(tif_condition_dir) if s.endswith('.tif'))
        except:
            print("\nCannot find tif file for %s" % condition)
            print(os.listdir(tif_condition_dir))

        fluor_names.append(os.path.join(fluordir, condition, tif_fn.strip('.tif') + '.CPfluor'))

    return fluor_names

def write_old_mapfile_for_processData(fluordir, mapfile, outfile):
    """
    Write a mapfile in the old format to be consistent with the requirement of `processData.py`
    Args:
        fluordir - str, absolute directory for CPfluor dirs, e.g. /.../data/fluor/
        mapfile - str
            a csv file with the column 'condition'
        outfile - str, name of the temporary output file
    """
    df = pd.read_csv(mapfile)
    conditions = df['condition'].tolist()
    
    lines = [fluordir, '.'] + conditions
    lines = [l + '\n' for l in lines]
    print(lines)
    
    with open(outfile, 'w+') as fh:
        fh.writelines(lines)

    print('Wrote to old style mapfile %s' % outfile)

def write_old_mapfile_from_fluordir(fluordir, outfile):
    """
    Write an old-fashioned mapfile with all the directories (conditions) in a clean fluordir.
    Args:
        fluordir - str
        outfile - str
    """
    conditions = os.listdir(fluordir)
    lines = [fluordir, '.'] + conditions
    lines = [l + '\n' for l in lines]
    
    with open(outfile, 'w+') as fh:
        fh.writelines(lines)
    
    print('Wrote to %s' % outfile)

def fill_sample_sheet(mapfile):
    """
    Fill tabular csv sample sheet.
    Args:
        mapfile - str
            a csv file with the column 'condition' 
            that's directory names e.g. '15_Green_15C'
    """
    df = pd.read_csv(mapfile)
    df['channel'] = df['condition'].apply(lambda s: s.split('_')[1])
    df['treatment'] = df['condition'].apply(lambda s: s.split('_')[2])
    df['is_temperature'] = df['treatment'].apply(lambda s: s.endswith('C'))

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

#fill_sample_sheet(r'config/nnnlib2b_map.csv')
#print(parse_mapfile(r'../config/test.map'))
#print(parse_fluorfiles_from_mapfile(r'../config/test.map'))
#datadir = r'/scratch/groups/wjg/kyx/NNNlib2b_Oct6/data/'
#fluor_names = get_fluor_names_from_mapfile(r'config/nnnlib2b_map.csv', fluordir=datadir+'fluor/', tifdir = datadir+'images/')
#print(fluor_names)
