import os
from scripts.util import *

####### SELECT CONFIG FILE HERE #######
configfile: "config/config_NNNlib2b_Oct6.yaml"
#######################################

# --- Define Global Variables --- #

datadir = config['datadir']
expdir = os.path.normpath(datadir + '/../') + '/'

# hardcoded tile numbers
TILES = ['tile%03d'%i for i in range(1,19)]
TILES_NO_ZERO_PAD = ['tile%d'%i for i in range(1,19)]

fluor_files = get_fluor_names_from_mapfile(config["mapfile"], config["tifdir"], config["fluordir"])

#wildcard_constraints:

# --- Define Required Output --- #

rule all:
    input: 
        config["seriesdir"]
        #expand(datadir + "filtered_tiles_libregion/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES),
        #fluor_files
        #datadir + "fluor/Green16_25/NNNlib2b_DNA_tile1_green_600ms_2011.10.22-16.51.13.953.CPfluor"
        #expand(expdir + "fig/fiducial/{tile}_Bottom_fiducial.png", tile=TILES)

#STRSTAMP, TILES = glob_wildcards(datadir + "tiles/{strstamp}_ALL_{tile}_Bottom.CPseq")

# --- Rules --- #

rule merge_fastqs_to_CPseq:
    input:
        r1 = config['fastq']['read1'],
        r2 = config['fastq']['read2']
    output:
        datadir + "sequence/ALL.CPseq"
    params:
        cluster_memory = "10G",
        cluster_time = "10:00:00"
    threads:
        2
    conda:
        "envs/ame.yml"
    shell:
        """
        python scripts/array_tools/CPscripts/mergeFastqReadsToCPseq.py -r1 {input.r1} -r2 {input.r2} -o {output} 
        """

rule split_CPseq:
    input:
        datadir + "sequence/ALL.CPseq"
    output:
        #directory(datadir + "tiles/")
        expand(datadir + "tiles/ALL_{tile}_Bottom.CPseq", tile=TILES)
    threads:
        1
    params:
        cluster_memory = "1G",
        cluster_time = "5:00:00",
        tiledir = datadir + "tiles/"
    conda:
        "envs/ame.yml"
    shell:
        """
        python scripts/array_tools/CPscripts/splitCPseqIntoTiles.py -o {params.tiledir} -s bottom {input}
        """

rule filter_tiles:
    input:
        expand(datadir + "tiles/ALL_{tile}_Bottom.CPseq", tile=TILES),
        config["FIDfilter"]
    output:
        expand(datadir + "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)
    params:
        tiledir = datadir + "tiles/",
        filteredtiledir = datadir + "filtered_tiles/",
        cluster_memory = "16G",
        cluster_time = "5:00:00"
    conda:
        "envs/ame.yml"
    #envmodules:
    #    "matlab"
    threads:
        8  
    shell:
        """
        module load matlab
        export MATLABPATH=/share/PI/wjg/lab/array_tools/CPscripts/:/share/PI/wjg/lab/array_tools/CPlibs/
        python scripts/array_tools/CPscripts/alignmentFilterMultiple.py -rd {params.tiledir} -f {config[FIDfilter]} -od {params.filteredtiledir} -gv /share/PI/wjg/lab/array_tools -n 18 
        """
rule filter_tiles_libregion:
    input:
        expand(datadir + "tiles/ALL_{tile}_Bottom.CPseq", tile=TILES),
        config["LibRegionFilter"]
    output:
        expand(datadir + "filtered_tiles_libregion/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)
    params:
        tiledir = datadir + "tiles/",
        filteredtiledir = datadir + "filtered_tiles_libregion/",
        cluster_memory = "16G",
        cluster_time = "5:00:00"
    conda:
        "envs/ame.yml"
    threads:
        8
    shell:
        """
        module load matlab
        export MATLABPATH=scripts/array_tools/CPscripts/:scripts/array_tools/CPlibs/
        python scripts/array_tools/CPscripts/alignmentFilterMultiple.py -rd {params.tiledir} -f {config[LibRegionFilter]} -od {params.filteredtiledir} -gv scripts/array_tools -n 18
        """

        
rule plot_fiducials:
    input:
        expand(datadir + "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)
    output:
        expand(expdir + "fig/fiducial/{tile}_Bottom_fiducial.png", tile=TILES)
    conda:
        "envs/plotting.yml"
    params:
        cluster_memory = "4G"
    threads:
        1
    script:
        "scripts/plot_seqs.py"

## quantify_images: quantify intensities in tif and write to CPfluor
## snakemake checks one tile per condition as input/output and submit one job per condition
rule quantify_images:
    input:
        image = config["tifdir"] + "{condition}/%s_{tile}_{channel}_600ms_{timestamp}.tif" % config["experimentName"],
        libregion = expand(datadir + "filtered_tiles_libregion/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)
    output:
        CPfluor = datadir + "fluor/{condition}/%s_{tile}_{channel}_600ms_{timestamp}.CPfluor" % config["experimentName"]#,
        #roff = expand(datadir + "roff/{condition}/%s_{tile}_{channel}_600ms_{timestamp}.roff")
    params:
        image_dir = config["tifdir"] + "{condition}/",
        seq_dir = datadir + "filtered_tiles_libregion/",
        fluor_dir = config["fluordir"] + "{condition}/",
        roff_dir = datadir + "roff/{condition}/",
        reg_subset = "LibRegion",
        log_dir = expdir + "log/quantify_image_{condition}.log",
        num_cores = "18",
        data_scaling = "MiSeq_to_TIRFStation1",
        cluster_memory = "40G",
        cluster_time = "15:00:00"
    threads:
        18
    conda:
        "envs/ame.yml"
    shell:
        """
        module load matlab
        export MATLABPATH=scripts/array_tools/CPscripts:scripts/array_tools/CPlibs
        python scripts/array_tools/CPscripts/quantifyTilesDownstream.py -id {params.image_dir} -ftd {params.seq_dir} -fd {params.fluor_dir} -rod {params.roff_dir} -n {params.num_cores} -rs {params.reg_subset} -sf {params.data_scaling} -gv scripts/array_tools/
        """


## combine_signal: Integrate and combine CPfluor files of different conditions into a single CPseries file per tile
rule combine_signal:
    input:
        fluor_files
    output:
        directory(config["seriesdir"])
    params:
        fluordir = config["fluordir"],
        libdata = config["libdata"],
        oldmapfile = datadir + 'tmp/' + config["experimentName"] + '.map',
        cluster_memory = "80G",
        cluster_time = "10:00:00",
        num_cores = "6"
    threads:
        6
    conda:
        "envs/py36.yml"
    shell:
        """
        python scripts/write_old_mapfile.py {params.fluordir} {config[mapfile]} {params.oldmapfile}
        python scripts/array_tools/bin_py3/processData.py -mf {params.oldmapfile} -od {output} --appendLibData {params.libdata} --num_cores {params.num_cores}
        """

"""

rule normalize:
    input:
        signal=""
        re=""
    shell:
        "python -m array_fitting_tools/bin/normalizeSeries -b {input.signal}.CPseries.gz -a {input.ref}.CPseries.gz"
"""