import pdb
import matplotlib.pyplot as plt
from tqdm import tqdm

def plot_species_per_tile(file_name, contains=None, color = 'b'):
    filename = file_name# + ".CPseq"

    # Using readlines()
    file1 = open(filename, 'r')
    Lines = file1.readlines()

    x_s = []
    y_s = []
    for i in tqdm(range(len(Lines))):
        line = Lines[i]
    # Look for lines that are the machine ID
        if line[:6] == 'M00653':
            line_split = line.split("\t") 
            if contains in line_split[1] or contains in line_split[2] or contains in line_split[4]:
                locations = line_split[0].split(":")
                x = int(locations[5])
                y = int(locations[6])*-1
                x_s.append(x)
                y_s.append(y)
    if color == "b":
        plt.scatter(x_s, y_s, s=1, c = color)
    else:
        plt.scatter(x_s, y_s, s=1, c = color)

def plot_tile(file_name):
    plt.figure(figsize=(14,11))
file_name = "/scratch/groups/wjg/mhinks/210612-miseq-8mers-allostry/filtered_tiles/JMHLB_ALL_tile004_Bottom_filtered"
#plot_tile(file_name, "GGGTCGGCTTCGCATATG", 'pink') # petcon3F
#plot_tile(file_name, "GGTGGATCAGGAGGTTCG", 'pink') # GS F
plot_tile(file_name, "AACAATTGCAGTGTT", 'fuchsia') # pho4
plot_tile(file_name, "CCTTTGT", 'fuchsia') #sox2sox2
#plot_tile(file_name, "GGGTCGGCTTCGCATATG", 'fuchsia') # petcon3 forward - will be all RNA except curtis's 2xFLAG
plot_tile(file_name, 'FID', 'b')# fiducial blue

plt.axis('off')
plt.savefig(file_name + "FID_sox2_.png", dpi=300, bbox_inches='tight')
