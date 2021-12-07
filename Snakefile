import os
from scripts.util import *

####### SELECT CONFIG FILE HERE #######
configfile: "config/config_NNNlib2b_20211022.yaml"
#######################################

# --- Define Global Variables --- #

datadir = config['datadir']
expdir = os.path.normpath(datadir + '/../') + '/'
sequencingResult = datadir + 'aligned/' + config["sequencingResult"]
normalizedSeries = datadir + "series_normalized/" + config["imagingExperiment"] + "_normalized.pkl"

# hardcoded tile numbers
TILES = ['tile%03d'%i for i in range(1,19)]
TILES_NO_ZERO_PAD = ['tile%d'%i for i in range(1,19)]

# == Decide which outputs to require base on progress of the experiment ==
assert config["processingType"] in ['pre-array', 'post-array']
if config["processingType"] == "pre-array":
    fluor_files = []
    requested_output = ["%s_STATS.csv" % sequencingResult.strip('.CPseq'),
                        expand(expdir + "fig/fiducial/{tile}_Bottom_fiducial.png", tile=TILES)]
elif config["processingType"] == "post-array":
    fluor_files = get_fluor_names_from_mapfile(config["mapfile"], config["tifdir"], config["fluordir"])
    requested_output = config["seriesdir"]#[fluor_files, config["seriesdir"]]
    
    if config["fitting"] == "NNN":
        requested_output.append(normalizedSeries)

#wildcard_constraints:

# --- Define Required Output --- #

rule all:
    input:
        # requested_output
        #"%s_STATS.csv" % sequencingResult.strip('.CPseq'),
        #config["sequencingResult"]#, #== Align sequencing data ==
        #expand(expdir + "fig/fiducial/{tile}_Bottom_fiducial.png", tile=TILES) #== Plot fiducials ==
        #expand(datadir + "filtered_tiles_libregion/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES), #== Filtered libregions ==
        #fluor_files
        #config["seriesdir"]
        #datadir + "fluor/Green16_25/NNNlib2b_DNA_tile1_green_600ms_2011.10.22-16.51.13.953.CPfluor" #== Example image quantification ==
        normalizedSeries


# --- Rules --- #

## unzip_fastq: in case the fastq files were not unzipped
rule unzip_fastq:
    input:
        r1 = config['fastq']['read1'] + ".gz",
        r2 = config['fastq']['read2'] + ".gz"
    output:
        r1 = config['fastq']['read1'],
        r2 = config['fastq']['read2']
    threads:
        1
    params:
        cluster_time = "02:00:00"
    shell:
        "gunzip {input.r1} {input.r2}"


## run_FLASH: align paired ends with FLASH
rule run_FLASH:
    input:
        r1 = config['fastq']['read1'],
        r2 = config['fastq']['read2']
    output:
        datadir + "FLASH_output/out.extendedFrags.fastq"
    params:
        outdir = datadir + "FLASH_output"
    threads:
        1
    shell:
        """
        cd {params[outdir]}
        {config[FLASHdir]}/FLASH-1.2.11/flash {input.r1} {input.r2}
        """    

rule convert_FLASH_to_CPseq:
    input:
        datadir + "FLASH_output/out.extendedFrags.fastq"
    output:
        datadir + "paired_reads/ConsensusPairedReads.CPseq"
    threads:
        1
    conda:
        "envs/align.yml"
    shell:
        "python scripts/convertFLASH_OutputToCPseq.py {input} {output}"



rule align_consensus_read_to_library:
    input:
        reads = datadir + "paired_reads/ConsensusPairedReads.CPseq",
        reference = config["referenceLibrary"],
        scoring_matrix = os.path.join(os.getcwd(), "data/reference/NUC.4.4") # need this to check existence of the matrix file
    output:
        sequencingResult
    threads:
        6
    params:
        fiveprimeRegion = config["refSeqContext"]["fiveprime"],
        threeprimeRegion = config["refSeqContext"]["threeprime"],
        cluster_memory = "90G",
        cluster_time = "24:00:00"    
    conda:
        "envs/align.yml"
    shell:
        """
        python scripts/matchConsensusReadsToLibrary.py --beam {config[alignBeam]} {input.reads} --library {input.reference} -o {output} --scoringMatrix {input.scoring_matrix} --fiveprimeRegion {params.fiveprimeRegion} --threeprimeRegion {params.threeprimeRegion}
        """

rule get_stats:
    input:
        sequencingResult
    output:
        "%s_STATS.csv" % sequencingResult.strip('.CPseq')
    threads:
        1
    conda:
        "envs/align.yml"
    shell:
        "python scripts/get_stats.py {input} {output}"


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
        "scripts/plotSeqs.py"


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


## write_old_mapfile: convert and prepare mapfile for the combine_signal step
rule write_old_mapfile:
    input:
        config['mapfile']
    output:
        oldmapfile = datadir + 'tmp/' + config["experimentName"] + '.map'
    params:
        fluordir = config["fluordir"],
        cluster_memory = "500M",
        cluster_time = "0:15:00"
    threads:
        1
    conda:
        "envs/py36.yml"
    shell:
        "python scripts/writeOldMapfile.py {params.fluordir} {config[mapfile]} {output.oldmapfile}"


## combine_signal: Integrate and combine CPfluor files of different conditions into a single CPseries file per tile
rule combine_signal:
    input:
        fluorfiles = fluor_files,
        oldmapfile = datadir + 'tmp/' + config["experimentName"] + '.map',
        libdata = config["sequencingResult"]
    output:
        directory(config["seriesdir"])
    params:
        cluster_memory = "80G",
        cluster_time = "10:00:00",
        num_cores = "6"
    threads:
        6
    conda:
        "envs/py36.yml"
    shell:
        """
        python scripts/array_tools/bin_py3/processData.py -mf {input.oldmapfile} -od {output} --appendLibData {input.libdata} --num_cores {params.num_cores}
        """

## normalize_signal: Given one merged CPseq file, normalize fluorescence signal for single cluster fit
rule normalize_signal:
    input:
        CPseries_file = config["seriesfile"],
        mapfile = config["mapfile"],
        annotation = config["referenceLibrary"]
    output:
        figdir = expdir + "fig/normalization",
        out_file = config["imagingExperiment"] + "_normalized.pkl",
        xdata_file = config["imagingExperiment"] + "_xdata.txt"
    params:
        green_norm_condition = "Green07_PostCy3",
        ext = ".pdf"
    threads:
        1
    script:
        "scripts/normalizeNNNlib2bSignal.py"

## fit_single_cluster
rule fit_single_cluster:
    input:
        normalized = normalizedSeries,
        xdata = config["imagingExperiment"] + "_xdata.txt"
    output:
        datadir + "fitted_single_cluster/" + config["imagingExperiment"] + ".CPfitted.gz"
    threads:
        12
    params:
        cluster_time = "20:00:00",
        cluster_memory = "32G"
    shell:
        "python scripts/nnn_fitting/singleClusterFits.py -b {input.normalized} -x {input.xdata} -o {output} --func melt_curve"