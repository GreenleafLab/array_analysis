import os

if __name__ == "__main__":
    flash_path = snakemake.config['FLASHdir']
    fastq_r1 = snakemake.config['fastq']['read1']
    fastq_r2 = snakemake.config['fastq']['read2']
    aligned_path = os.path.join(snakemake.config['datadir'], 'aligned')

    os

