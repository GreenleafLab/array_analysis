import sys
from util import *

if __name__ == "__main__":
    fluordir = sys.argv[1]
    mapfile = sys.argv[2]
    outdir = sys.argv[3]
    write_old_mapfile_for_processData(fluordir, mapfile, outdir)
