import os
# import scripts.getSnakeConfig as snakeconfig

# configfile: "config.yaml"

# fluorfile,  = snakeconfig.parse_mapfile(config['mapfile'])
# CHANNELS = ["Green", "Red"]

rule CPfluor2CPseries:
    input:
        "data/CPfluor/test.CPfluor"
    output:
        "data/CPseries/test.CPseries"
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