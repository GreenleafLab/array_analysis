import os
from scripts.util import *

####### SELECT CONFIG FILE HERE #######
configfile: "config/config_rf009.yaml"
#######################################

# --- Define Global Variables --- #

datadir = config['datadir']
expdir = os.path.normpath(os.path.join(datadir, '..'))
sequencingResult = config["sequencingResult"]
normalizedSeries = os.path.join(datadir, "series_normalized/",  config["imagingExperiment"] + "_normalized.pkl")
fittedVariant = os.path.join(datadir, "fitted_variant/",  config["imagingExperiment"] + ".CPvariant.gz")
# hardcoded tile numbers
TILES = ['tile%03d'%i for i in range(1,19)]
TILES_NO_ZERO_PAD = ['tile%d'%i for i in range(1,19)]

# == Decide which outputs to require base on progress of the experiment ==
assert config["processingType"] in ['pre-array', 'post-array']
if config["processingType"] == "pre-array":
    fluor_files = []
    requested_output = ["%s_STATS.csv" % sequencingResult.replace('.CPseq', ''),
                        expand(os.path.join(expdir, "fig/fiducial/{tile}_Bottom_fiducial.png"), tile=TILES)]
elif config["processingType"] == "post-array":
    fluor_files = get_fluor_names_from_mapfile(config["mapfile"], config["tifdir"], config["fluordir"])
    if config["fitting"] == "NNN":
        requested_output = fittedVariant
    else:
        requested_output = config["seriesfile"]
    
container: "library://kyx/array_analysis/mamba:latest"
#wildcard_constraints:

# --- Define Required Output --- #

rule all:
    input:
        requested_output


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
        os.path.join(datadir, "FLASH_output/out.extendedFrags.fastq")
    params:
        outdir = os.path.join(datadir, "FLASH_output")
    threads:
        1
    shell:
        """
        cd {params[outdir]}
        {config[FLASHdir]}/FLASH-1.2.11/flash {input.r1} {input.r2}
        """    

rule convert_FLASH_to_CPseq:
    input:
        os.path.join(datadir, "FLASH_output/out.extendedFrags.fastq")
    output:
        os.path.join(datadir, "paired_reads/ConsensusPairedReads.CPseq")
    params:
        cluster_memory = "16G"
    threads:
        1
    conda:
        "envs/align.yml"
    shell:
        "python3 scripts/convertFLASH_OutputToCPseq.py {input} {output}"



rule align_consensus_read_to_library:
    input:
        reads = os.path.join(datadir, "paired_reads/ConsensusPairedReads.CPseq"),
        reference = config["referenceLibrary"],
        scoring_matrix = os.path.join(os.getcwd(), "data/reference/NUC.4.4") # need this to check existence of the matrix file
    output:
        seq_result = sequencingResult,
        cluster_annot = sequencingResult.replace('.CPseq', '.CPannot')
    threads:
        6
    params:
        fiveprimeRegion = config["refSeqContext"]["fiveprime"],
        threeprimeRegion = config["refSeqContext"]["threeprime"],
        variant_col = config['variantCol'],
        cluster_memory = "90G",
        cluster_time = "24:00:00"    
    conda:
        "envs/align.yml"
    shell:
        """
        python3 scripts/matchConsensusReadsToLibrary.py --beam {config[alignBeam]} {input.reads} --library {input.reference} -o {output.seq_result} --clusterAnnotation {output.cluster_annot} --scoringMatrix {input.scoring_matrix} --fiveprimeRegion {params.fiveprimeRegion} --threeprimeRegion {params.threeprimeRegion}
        """

rule get_stats:
    input:
        sequencingResult
    output:
        "%s_STATS.csv" % sequencingResult.strip('.CPseq')
    params:
        cluster_memory = "16G"
    threads:
        1
    conda:
        "envs/align.yml"
    shell:
        "python3 scripts/get_stats.py {input} {output}"


rule merge_fastqs_to_CPseq:
    input:
        r1 = config['fastq']['read1'],
        r2 = config['fastq']['read2']
    output:
        os.path.join(datadir, "sequence/ALL.CPseq")
    params:
        cluster_memory = "10G",
        cluster_time = "10:00:00"
    threads:
        2
    conda:
        "envs/py36.yml"
    shell:
        """
        python3 array_tools/CPscripts/mergeFastqReadsToCPseq.py -r1 {input.r1} -r2 {input.r2} -o {output} 
        """

rule split_CPseq:
    input:
        os.path.join(datadir, "sequence/ALL.CPseq")
    output:
        expand(os.path.join(datadir, "tiles/ALL_{tile}_Bottom.CPseq"), tile=TILES)
    threads:
        1
    params:
        cluster_memory = "1G",
        cluster_time = "5:00:00",
        tiledir = os.path.join(datadir, "tiles/")
    conda:
        "envs/py36.yml"
    shell:
        """
        python3 array_tools/CPscripts/splitCPseqIntoTiles.py -o {params.tiledir} -s bottom {input}
        """

rule filter_tiles:
    input:
        expand(os.path.join(datadir, "tiles/ALL_{tile}_Bottom.CPseq"), tile=TILES),
        config["FIDfilter"]
    output:
        expand(os.path.join(datadir, "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq"), tile=TILES)
    params:
        tiledir = os.path.join(datadir, "tiles/"),
        filteredtiledir = os.path.join(datadir, "filtered_tiles/"),
        cluster_memory = "24G",
        cluster_time = "5:00:00"
    conda:
        "envs/py36.yml"
    #envmodules:
    #    "matlab"
    threads:
        8  
    shell:
        """
        module load matlab
        export MATLABPATH=array_tools/CPscripts/:array_tools/CPlibs/
        python3 array_tools/CPscripts/alignmentFilterMultiple.py -rd {params.tiledir} -f {config[FIDfilter]} -od {params.filteredtiledir} -gv array_tools -n 18 
       """
       
rule filter_tiles_libregion:
    input:
        expand(os.path.join(datadir, "tiles/ALL_{tile}_Bottom.CPseq"), tile=TILES),
        config["LibRegionFilter"]
    output:
        expand(os.path.join(datadir, "filtered_tiles_libregion/ALL_{tile}_Bottom_filtered.CPseq"), tile=TILES)
    params:
        tiledir = os.path.join(datadir, "tiles/"),
        filteredtiledir = os.path.join(datadir, "filtered_tiles_libregion/"),
        cluster_memory = "32G",
        cluster_time = "5:00:00"
    conda:
        "envs/py36.yml"
    threads:
        8
    shell:
        """
        module load matlab
        export MATLABPATH=array_tools/CPscripts/:array_tools/CPlibs/
        python3 array_tools/CPscripts/alignmentFilterMultiple.py -rd {params.tiledir} -f {config[LibRegionFilter]} -od {params.filteredtiledir} -gv array_tools -n 18
        """

        
rule plot_fiducials:
    input:
        expand(os.path.join(datadir, "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq"), tile=TILES)
    output:
        expand(os.path.join(expdir, "fig/fiducial/{tile}_Bottom_fiducial.png"), tile=TILES)
    conda:
        "envs/plotting.yml"
    params:
        cluster_memory = "4G"
    threads:
        1
    script:
        "scripts/plotSeqs.py"


# ## quantify_images: quantify intensities in tif and write to CPfluor
# ## snakemake checks one tile per condition as input/output and submit one job per condition
# rule quantify_images:
#     input:
#         image = os.path.join(config["tifdir"], "{condition}/%s_{tile}_{channel}_{exposure}_{timestamp}.tif") % config["prefix"],
#         libregion = expand(os.path.join(datadir, "filtered_tiles_libregion/ALL_{tile}_Bottom_filtered.CPseq"), tile=TILES)
#     output:
#         CPfluor = os.path.join(config["fluordir"], "{condition}/%s_{tile}_{channel}_{exposure}_{timestamp}.CPfluor") % config["prefix"]#,
#     params:
#         image_dir = os.path.join(config["tifdir"],  "{condition}/"),
#         seq_dir = os.path.join(datadir, "filtered_tiles_libregion/"),
#         fluor_dir = os.path.join(config["fluordir"],  "{condition}/"),
#         roff_dir = os.path.join(datadir, "roff/{condition}/"),
#         reg_subset = "LibRegion",
#         log_dir = os.path.join(expdir,  "log/quantify_image_{condition}.log"),
#         num_cores = "18",
#         data_scaling = "MiSeq_to_TIRFStation1",
#         cluster_memory = "40G",
#         cluster_time = "15:00:00"
#     threads:
#         18
#     conda:
#         "envs/py36.yml"
#     shell:
#         """
#         module load matlab
#         export MATLABPATH=array_tools/CPscripts:array_tools/CPlibs
#         python3 array_tools/CPscripts/quantifyTilesDownstream.py -id {params.image_dir} -ftd {params.seq_dir} -fd {params.fluor_dir} -rod {params.roff_dir} -n {params.num_cores} -rs {params.reg_subset} -sf {params.data_scaling} -gv array_tools/
#         """


# ## write_old_mapfile: convert and prepare mapfile for the combine_signal step
# rule write_old_mapfile:
#     input:
#         config['mapfile']
#     output:
#         oldmapfile = os.path.join(datadir, 'tmp/', config["imagingExperiment"] + '.map')
#     params:
#         fluordir = config["fluordir"],
#         cluster_memory = "500M",
#         cluster_time = "0:15:00"
#     threads:
#         1
#     conda:
#         "envs/py36.yml"
#     shell:
#         "python3 scripts/writeOldMapfile.py {params.fluordir} {config[mapfile]} {output.oldmapfile}"


# ## combine_signal: Integrate and combine CPfluor files of different conditions into a single CPseries file per tile
# rule combine_signal:
#     input:
#         fluorfiles = fluor_files,
#         oldmapfile = os.path.join(datadir, 'tmp/', config["imagingExperiment"] + '.map'),
#         libdata = sequencingResult
#     output:
#         get_series_tile_filenames(config["seriesdir"], config["prefix"])
#     params:
#         output_directory = config["seriesdir"],
#         cluster_memory = "80G",
#         cluster_time = "00:30:00",
#         num_cores = "3"
#     threads:
#         3
#     conda:
#         "envs/py36.yml"
#     shell:
#         """
#         python3 array_tools/bin_py3/processData.py -mf {input.oldmapfile} -od {params.output_directory} --appendLibData {input.libdata} --num_cores {params.num_cores}
#         """

# ## concat_tiles_signal: Concatenate one CPseries file per tile into one big CPseries file with all tiles
# rule concat_tiles_signal:
#     input:
#         series = ancient(get_series_tile_filenames(config["seriesdir"], config["prefix"])),
#         mapfile = ancient(config["mapfile"])
#     output:
#         config["seriesfile"]
#     conda:
#         "envs/py36.yml"
#     shell:
#         """
#         python3 scripts/concatTilesSignal.py --tiles {input.series} -o {output} -m {input.mapfile}
#         """

## normalize_signal: Given one merged CPseq file, normalize fluorescence signal for single cluster fit
rule normalize_signal:
    input:
        CPseries_file = ancient(config["seriesfile"]),
        mapfile = ancient(config["mapfile"]),
        annotation = ancient(config["referenceLibrary"])
    output:
        out_file = os.path.join(datadir, 'series_normalized/',  "%s_normalized.pkl" % config["imagingExperiment"]),
        xdata_file = os.path.join(datadir, 'series_normalized/', "%s_xdata.txt" % config["imagingExperiment"])
    params:
        figdir = os.path.join(expdir,  "fig/normalization_%s/"%config["imagingExperiment"]),
        green_norm_condition = config["greenNormCondition"],
        ext = ".pdf",
        variant_col = config["variantCol"],
        smooth = config['smooth'] # savgol_{window_length}_{polyorder} e.g. savgol_7_2 or 'None'
    conda:
        "envs/fitting.yml"    
    threads:
        1
    script:
        "scripts/normalizeNNNlib2bSignal.py"

## fit_single_cluster
rule fit_single_cluster:
    input:
        normalized = normalizedSeries,
        xdata = os.path.join(datadir, 'series_normalized/', "%s_xdata.txt" % config["imagingExperiment"]),
        mapfile = config["mapfile"]
    output:
        os.path.join(datadir, "fitted_single_cluster/", "%s.CPfitted.gz" % config["imagingExperiment"])
    threads:
        18
    params:
        cluster_time = "48:00:00",
        cluster_memory = "32G"
    conda:
        "envs/fitting.yml"
    shell:
        "python3 scripts/nnn_fitting/singleClusterFits.py --parallel -b {input.normalized} -x {input.xdata} -o {output} --mapfile {input.mapfile}"


## bootstrap_variant_median
rule bootstrap_variant_median:
    input:
        cf = os.path.join(datadir, "fitted_single_cluster/",  config["imagingExperiment"] + ".CPfitted.gz"),
        annotation = sequencingResult.replace('.CPseq', '.CPannot')
    output:
        variant = os.path.join(datadir, "fitted_single_cluster/",  config["imagingExperiment"] + ".CPvariant"),
        good_clusters = os.path.join(datadir, "fitted_single_cluster/",  config["imagingExperiment"] + "_good_cluster_ind.txt")
    params:
        p = "dH Tm",
        n_samples = "100",
        good_fit = config["query"]["singleCluster"].replace(' ', ''),
        vc = config["variantCol"],
        cluster_time = "02:00:00",
        cluster_memory = "8G"
    threads:
        12
    conda:
        "envs/fitting.yml"
    shell:
        """
        python3 scripts/nnn_fitting/bootStrapFitFile.py -cf {input.cf} -a {input.annotation} -g {output.good_clusters}\
                -p {params.p} -vc {params.vc} --query {params.good_fit} --n_samples {params.n_samples}
        """


## fit_fmax_fmin_distribution
rule fit_fmax_fmin_distribution:
    input:
        vf = os.path.join(datadir, "fitted_single_cluster/",  config["imagingExperiment"] + ".CPvariant")
    output:
        fm = os.path.join(datadir, "fitted_fmax_fmin/%s-fmax_fmin.json" % config["imagingExperiment"]),
        plots = directory(os.path.join(expdir,  "fig/fmax_fmin_%s/"%config["imagingExperiment"])),
        # plots = expand(expdir,  "fig/fmax_fmin/{plotname}"%config["experimentName"], plotname=['fmax_vs_dG_init.pdf', 'fmin_vs_dG_init.pdf'])
    params:
        figdir = os.path.join(expdir,  "fig/fmax_fmin_%s/"%config["imagingExperiment"]),
        fmax_q = config["query"]["fmaxVariant"].replace(' ', ''),
        fmin_q = config["query"]["fminVariant"].replace(' ', ''),
        variant_q = config["query"]["variant"].replace(' ', '')
    threads:
        1
    conda:
        "envs/fitting.yml"
    shell:
        """
        python3 scripts/nnn_fitting/findFmaxFminDist.py -vf {input.vf} -o {output.fm} --figdir {params.figdir}\
                -fmaxq {params.fmax_q} -fminq {params.fmin_q} --variant_filter {params.variant_q}
        """


## fit_refine_variant: refine fit at the variant level using estimated fmax and fmin for those not reaching them
rule fit_refine_variant:
    input:
        cluster = normalizedSeries,
        variant = os.path.join(datadir, "fitted_single_cluster/",  config["imagingExperiment"] + ".CPvariant"),
        xdata = os.path.join(datadir, "series_normalized/",  config["imagingExperiment"] + "_xdata.txt"),
        mapfile = config["mapfile"],
        fm = os.path.join(datadir, "fitted_fmax_fmin/%s-fmax_fmin.json" % config["imagingExperiment"]),
        annotation = sequencingResult.replace('.CPseq', '.CPannot')
    output:
        fitted = os.path.join(datadir, "fitted_variant/",  config["imagingExperiment"] + ".CPvariant.gz")
    params:
        figdir = os.path.join(expdir,  "fig/%s-fit_refine_variant/"%config["imagingExperiment"]),
        p = "dH Tm",
        n_bootstraps = "100",
        vc = config["variantCol"],
        good_clusters = os.path.join(datadir, "fitted_single_cluster/", config["imagingExperiment"] + "_good_cluster_ind.txt"),
        variant_q = config["query"]["variant"].replace(" ", ""),
        cluster_time = "48:00:00",
        cluster_memory = "48G"
    threads:
        20
    conda:
        "envs/fitting.yml"
    shell:
        """
        python3 scripts/nnn_fitting/refineVariantFits.py --parallel -b {input.cluster} -vf {input.variant} -x {input.xdata}\
        --mapfile {input.mapfile} --cluster_annot {input.annotation} --fmax_fmin {input.fm} -o {output.fitted} --figdir {params.figdir}\
        --param {params.p} --variant_col {params.vc} --n_bootstraps {params.n_bootstraps}
        """