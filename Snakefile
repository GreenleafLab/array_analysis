import scripts.getSnakeConfig as snakeconfig

configfile: "config.yaml"

fluordir, CONDITIONS = snakeconfig.parse_mapfile(config['mapfile'])
CHANNELS = ["Green", "Red"]

rule CPfluor2CPseries:
    input:
        

rule normalize:
    input:
        signal=""
        re=""
    shell:
        "python -m array_fitting_tools/bin/normalizeSeries -b {input.signal}.CPseries.gz -a {input.ref}.CPseries.gz"