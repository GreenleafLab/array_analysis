import os
import scripts.getSnakeConfig as snakeconfig

#configfile: "config/config.yaml"

fluorfiles, seriesfiles = snakeconfig.parse_fluorfiles_from_mapfile('config/test.map')

#wildcard_constraints:


rule all:
    input: seriesfiles

rule CPfluor2CPseries:
    input:
        "{basedir}/fluor_dir/{round}/{tile}.CPfluor"
    output:
        "{basedir}/series_dir/{round}/{tile}.CPseries"
    script:
        "scripts/CPfluor2CPseries.py"

"""
rule normalize:
    input:
        signal=""
        re=""
    shell:
        "python -m array_fitting_tools/bin/normalizeSeries -b {input.signal}.CPseries.gz -a {input.ref}.CPseries.gz"
"""