import os

def parse_mapfile(mapfile):
    """
    Args:
        mapfile - str
    Returns:
        fluordir - str
        condition - List[str]
    """
    with open(mapfile, 'r') as mf:
        mapinfo = mf.readlines()
    
    fluordir = mapinfo[0]
    # channel = [s.split("_")[1] for s in mapinfo[1:]]
    condition = [s.split("_")[-1].strip('\n') for s in mapinfo[1:]]

    return fluordir, condition

fluordir, condition = parse_mapfile('../config/nnnlib2.map')