import matplotlib.pyplot as plt
from tqdm import tqdm

def plot_species_per_tile(file_name, contains=None, color = 'b'):
    """
    Args:
        file_name - str, name of 1 filtered tiles CPseq file
    """
    # Using readlines()
    file1 = open(file_name, 'r')
    Lines = file1.readlines()

    x_s = []
    y_s = []
    for i in tqdm(range(len(Lines))):
        line = Lines[i]
        # Look for lines that are the matching ID
        if line[:6] == 'M03796':
            line_split = line.split("\t") 
            if contains in line_split[1] or contains in line_split[2] or contains in line_split[4]:
                locations = line_split[0].split(":")
                x = int(locations[5])
                y = int(locations[6])*-1
                x_s.append(x)
                y_s.append(y)
    
    
    print(x_s)
    print('here')
    plt.scatter(x_s, y_s, s=1, c = color)


def plot_fiducial_in_tile(CPseq_name, fig_name):
    """
    Args:
        CPseq_name - str, name of 1 filtered tiles CPseq file
        fig_name - str, output
    """
    plt.figure(figsize=(14,11))

    # plot one or multiple species here
    plot_species_per_tile(CPseq_name, 'FID', 'k')# fiducial blue

    plt.axis('off')
    plt.savefig(fig_name, dpi=300, bbox_inches='tight')


def plot_fiducial_all_tiles(CPseq_names, fig_names):
    """
    Args:
        CPseq_names - List[str]
        fig_names - List[str]
    """
    for seqfile,figfile in zip(CPseq_names, fig_names):
        #fig_name = os.path.join(fig_dir, os.path.split(seqfile)[1].replace('.CPseq', '.png'))
        plot_fiducial_in_tile(seqfile, figfile)

if __name__ == "__main__":
    plot_fiducial_all_tiles(snakemake.input[0], snakemake.output[0])
