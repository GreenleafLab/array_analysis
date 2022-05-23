import argparse
import numpy as np
import pandas as pd
import os,sys

def write_to_hdf5(series_tiles, h5nm, mapfile, clean=True):
        
    conditions = pd.read_csv(mapfile)['condition'].tolist()

    with pd.HDFStore(h5nm) as h5:
        for i in range(1, 19):
            print('Reading tile %03d' % i)
            tile = pd.read_csv(series_tiles[i - 1])
            
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


def write_to_pickle(series_tiles, pklnm, mapfile, clean=True):
    
    conditions = pd.read_csv(mapfile)['condition'].tolist()

    tiles = []
    for i in range(1,19):
        print('Reading tile %03d' % i)
        tile = pd.read_csv(series_tiles[i - 1])

        if all(condition in tile.columns for condition in conditions):
        # Check if all conditions are in the tile CPseries file

            if clean:
                tile.dropna(how='all', subset=conditions, inplace=True)
                tile.drop_duplicates(subset='clusterID', keep=False, ignore_index=True, inplace=True)

            tiles.append(tile)
        else:
            print('Tile CPseries file %s dropped for missing condition(s)'%series_tiles[i - 1])

    df = pd.concat(tiles, axis=0)
    df.to_pickle(pklnm)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='combine signal data from all tiles into one tile')    
    parser.add_argument('-t', '--tiles', nargs='+',
                    help='list of CPseries files, each for a tile')
    parser.add_argument('-o', '--output', help='name of the output file, could be either .h5 or .pkl')
    parser.add_argument('-m', '--mapfile', help='a csv file with the col `condition` for filtering nans')
    
    args = parser.parse_args()
    
    if args.output.endswith('.h5'):
        write_to_hdf5(series_tiles=args.tiles, pklnm=args.output, mapfile=args.mapfile, clean=True)
    elif args.output.endswith('.pkl'):
        write_to_pickle(series_tiles=args.tiles, pklnm=args.output, mapfile=args.mapfile, clean=True)
    else:
        print('Output file must be .h5 or .pkl')
