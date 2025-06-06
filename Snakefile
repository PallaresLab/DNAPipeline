import os
shell.executable("bash")

from snakemake.utils import min_version
min_version("8.10.6")

configfile: "config.yaml"
log_dir = config['output_dir'] + "/logs"
working_dir = config['output_dir'] + "/working"
ref_genome_name=os.path.basename(config['reference_genome'])
sample_files = snakemake.utils.listfiles(config["fastq_dir"]+"/{sample}.fastq.gz")
samples = dict((y[0], x) for x, y in sample_files)

ALL_SAMPLES = glob_wildcards(config['fastq_dir']+"/{sample}.fastq.gz").sample
PE_SAMPLES = glob_wildcards(config['fastq_dir']+"/{sample}_R1_001.fastq.gz").sample

def get_paired_fastq(wildcards):  
    paired_files = [samples[i] for i in samples.keys()
                if ("R1" in i or "R2" in i) and i.startswith(f"{wildcards.sample}_")]
    if len(paired_files) == 2 :
        return sorted(paired_files)
    else:
        raise ValueError(f"Error in matched pairs {wildcards.sample}")
        


    
rule all:
    input:
        config['output_dir']+"/multiqc/multiqc.html",
        config['output_dir']+"/multiqc/QC_table.csv"
        

    

rule fastqc:
    input:
        config['fastq_dir']+"/{sample}.fastq.gz",
        
    output:
        zip = working_dir+"/fastqc/{sample}_fastqc.zip",
        html = temp(working_dir+"/fastqc/{sample}_fastqc.html")
       
    log:
        log_dir + "/fastqc/{sample}.log",

    params:
        outdir = working_dir+"/fastqc",
    
        
    threads: 32
       
    conda:
        "envs/fastqc.yml"
      
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --threads {threads} {input} --outdir {params.outdir} >{log} 2>&1 "

       

rule trim:
    input:
        get_paired_fastq
       
    output:
        out_dir=temp(directory(working_dir+"/trimmed/{sample}")),
        R1 = temp(working_dir+"/trimmed/{sample}/{sample}_R1_001_val_1.fq.gz"),
        R2 = temp(working_dir+"/trimmed/{sample}/{sample}_R2_001_val_2.fq.gz"),
        log=log_dir+"/trimmed/{sample}.log"

    params:
        outdir=working_dir+"/trimmed/{sample}",
       
    threads: 32
    
    conda:
        "envs/trim_galore.yml"
       
    shell:
        "mkdir -p {params.outdir} && "
        "trim_galore --paired {input} --fastqc --output_dir {params.outdir} --cores {threads} "
        ">{output.log} 2>&1 "



        


rule bwa_index:
    input:
        genome = config['reference_genome']

    output:
        file = working_dir+"/bwa_index/"+ref_genome_name
        
    log:
        log_dir+"/bwa_index/index.log"
        
    params:
        outdir = working_dir+"/bwa_index"
        
    conda:
        "envs/bwa.yml"
        
    shell:
        "mkdir -p {params.outdir} &&"
        "cp {input.genome} {params.outdir} &&"
        "bwa-mem2 index {output.file} "
        ">{log} 2>&1"


rule bwa:
    input:
        R1 = working_dir+"/trimmed/{sample}/{sample}_R1_001_val_1.fq.gz",
        R2 = working_dir+"/trimmed/{sample}/{sample}_R2_001_val_2.fq.gz",
        genome=working_dir+"/bwa_index/"+ref_genome_name
        
    output:
        file = temp(working_dir+"/bwa/{sample}.bam"),
        stats = working_dir+"/bwa/{sample}.txt"
        
    log:
        log_dir+"/bwa/{sample}.log"
        
    params:
        outdir=working_dir+"/bwa"
        
    threads: 32
        
    conda:
        "envs/bwa.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "bwa-mem2 mem -M -t {threads} -R $(bash get_RG.sh {input.R1}) {input.genome} {input.R1} {input.R2} | samtools sort -@{threads} | samtools view -b -F 256 -f 2 -o {output.file} "# keep only paired mapped reads
        ">{log} 2>&1 && "
        "samtools stats {output.file} >{output.stats}"


        
rule picard:
    input:
        working_dir+"/bwa/{sample}.bam",

    output:
        file=working_dir+"/picard/{sample}_deduplicated.bam",
        metrics=working_dir+"/picard/{sample}_deduplicated.metrics.txt"
        
    log:
        log_dir+"/picard/{sample}.log"
        
    params:
        outdir=working_dir+"/picard"
        
    conda:
        "envs/picard.yml"  

    shell:
        "mkdir -p {params.outdir} &&"
        "picard MarkDuplicates REMOVE_DUPLICATES=true I={input} O={output.file} M={output.metrics} CREATE_INDEX=true"
        ">{log} 2>&1"
   
    
    
rule multiqc:
    input: 
        expand(working_dir+"/fastqc/{sample}_fastqc.zip",
               sample=ALL_SAMPLES),
        expand(working_dir+"/trimmed/{sample}",
               sample=PE_SAMPLES),
        expand(log_dir+"/trimmed/{sample}.log",
               sample=PE_SAMPLES),
        expand(working_dir+"/bwa/{sample}.txt",
               sample=PE_SAMPLES),    
        expand(working_dir+"/picard/{sample}_deduplicated.metrics.txt",
               sample=PE_SAMPLES),
               
          
    output:
        html=config['output_dir']+"/multiqc/multiqc.html",
        folder=directory(config['output_dir']+"/multiqc/multiqc_data")
        
    log:
        log_dir+"/multiqc/multiqc.log"
        
    params:
        outdir=config['output_dir']+"/multiqc"

    conda:
        "envs/multiqc.yml"
        
    shell:
        "mkdir -p {params.outdir} && "
        "multiqc --config multiqc.yaml -o {params.outdir} -n multiqc.html {input} "
        ">{log} 2>&1"


rule QC_table:
    input:
        config['output_dir']+"/multiqc/multiqc_data"
    output:
        config['output_dir']+"/multiqc/QC_table.csv"
    shell:
        "python DNA_QCtable.py -m {input} -o {output} "
