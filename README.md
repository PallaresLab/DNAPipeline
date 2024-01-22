# Paired-end DNA pipeline

## Description
This snakemake pipeline is designed for paired-end NGS DNA

### Inputs

*   Path of fastqc files
*   Name of outout folder 
*   Reference genome 


### Outputs

*   `multiqc_data\` - Dictionary containing the summary results of all the tools, inculde multiqc.html
*   `logs\` - Directory of log files for each job, check here first if you run into errors
*   `working\` - Directory containing intermediate files for each job

### Workflow

1.  **QC--fastqc 
2.  **Trimming--trim galore
3.  **QC--fastqc
4.  **Align--bwa
5.  **Sort--samtools
6.  **Deduplicate--picard
7.  **Summary--multiqc


## Setup environment

1.  Install conda

2.  Create a new enviroment

    ```bash
    conda create -n <project_name>
    ```

3.  Activate the environment

    ```bash
    conda activate <project_name>
    ```
    
4.  Enable the [Bioconda](https://bioconda.github.io/#using-bioconda) channel

    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

5. Install snakemake

    ```bash
    conda install snakemake
    ```

## Run workflow

1.  Clone workflow into working directory

    ```bash
    git clone <repo> <dir>
    cd <dir>
    ```

2.  Edit configuration files
    change the path of fastq_dir, output_dir, reference_genome in "config.yaml"


3.  Execute the workflow

    ```bash
    snakemake --configfile "config.yaml" --use-conda  --cores N
    ```

