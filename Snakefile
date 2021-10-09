import os
#import scripts.getSnakeConfig as snakeconfig

####### SELECT CONFIG FILE HERE #######
configfile: "config/config_NNNlib2b_Oct6.yaml"
#######################################

DATADIR = config['datadir']
# hardcoded tile numbers
TILES = ['tile%03d'%i for i in range(1,20)]

#_,_,_,ROUNDS = snakeconfig.parse_mapfile('config/nnnlib2.map')
#fluorfiles, seriesfiles = snakeconfig.parse_fluorfiles_from_mapfile('config/nnnlib2.map')

#wildcard_constraints:

rule all:
    input: expand(DATADIR + "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)


#STRSTAMP, TILES = glob_wildcards(DATADIR + "tiles/{strstamp}_ALL_{tile}_Bottom.CPseq")

rule merge_fastqs_to_CPseq:
    input:
        r1 = config['fastq']['read1'],
        r2 = config['fastq']['read2']
    output:
        DATADIR + "sequence/ALL.CPseq"
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
        DATADIR + "sequence/ALL.CPseq"
    output:
        #directory(DATADIR + "tiles/")
        expand(DATADIR + "tiles/ALL_{tile}_Bottom.CPseq", tile=TILES)
    threads:
        1
    params:
        cluster_memory = "1G",
        cluster_time = "5:00:00",
        tiledir = DATADIR + "tiles/"
    conda:
        "envs/ame.yml"
    shell:
        """
        python scripts/array_tools/CPscripts/splitCPseqIntoTiles.py -o {params.tiledir} -s bottom {input}
        """

rule filter_tiles:
    input:
        expand(DATADIR + "tiles/ALL_{tile}_Bottom.CPseq", tile=TILES)
    output:
        expand(DATADIR + "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)
    params:
        tiledir = DATADIR + "tiles/",
        filteredtiledir = DATADIR + "filtered_tiles/",
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

rule plot_fiducials:
    input:
        expand(DATADIR + "filtered_tiles/ALL_{tile}_Bottom_filtered.CPseq", tile=TILES)
    output:
        expand(DATADIR + "fig/{tile}_Bottom_fiducial.png", tile=TILES)
    conda:
        "envs/plotting.yml"
    threads:
        1
    script:
        "scripts/plot_seqs.py"

"""
rule integrate_signal:
    input:
        "{datadir}/fluor/{round}/{tile}.CPfluor"
    output:
        "{datadir}/signal/{round}/{tile}.CPseries"
    script:
        "scripts/integrateSignal.py"

rule merge_signal_2_series:
    input:
        "{datadir}/signal/{round}/"

rule normalize:
    input:
        signal=""
        re=""
    shell:
        "python -m array_fitting_tools/bin/normalizeSeries -b {input.signal}.CPseries.gz -a {input.ref}.CPseries.gz"
"""