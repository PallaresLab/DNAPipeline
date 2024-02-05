import os
configfile: "config.yaml"
log_dir = config['output_dir'] + "/logs"
working_dir = config['output_dir'] + "/working"
ref_genome_name=os.path.basename(config['reference_genome'])
SAMPLES=glob_wildcards(config['fastq_dir']+"/{sample}_R1_001.fastq.gz").sample

rule all:
    input:
        config['output_dir']+"/multiqc/multiqc.html"
        


rule fastqc:
    input:
        R1= config['fastq_dir']+"/{sample}_R1_001.fastq.gz",
        R2= config['fastq_dir']+"/{sample}_R2_001.fastq.gz"
        
    output:
        working_dir+"/fastqc/{sample}_R1_001_fastqc.zip",
        working_dir+"/fastqc/{sample}_R2_001_fastqc.zip"
       
    log:
        log1=log_dir+"/fastqc/{sample}_R1_001.log",
        log2=log_dir+"/fastqc/{sample}_R2_001.log"
      
    params:
        outdir=working_dir+"/fastqc"
        
    threads: 32
       
    conda:
        "envs/fastqc.yml"
      
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --threads {threads} {input.R1} --outdir {params.outdir} >{log.log1} 2>&1 &&"
        "fastqc --threads {threads} {input.R2} --outdir {params.outdir} >{log.log2} 2>&1 "
       

rule trim:
    input:
        R1=config['fastq_dir']+"/{sample}_R1_001.fastq.gz",
        R2=config['fastq_dir']+"/{sample}_R2_001.fastq.gz"

    output:
        working_dir+"/trimmed/{sample}_R1_001_val_1.fq.gz",
        working_dir+"/trimmed/{sample}_R2_001_val_2.fq.gz",
        working_dir+"/trimmed/{sample}_R1_001.fastq.gz_trimming_report.txt",
        working_dir+"/trimmed/{sample}_R2_001.fastq.gz_trimming_report.txt",
        working_dir+"/trimmed/{sample}_R1_001_val_1_fastqc.zip",
        working_dir+"/trimmed/{sample}_R2_001_val_2_fastqc.zip"
        
    log:
       log_dir+"/trimmed/{sample}.log"
       
    params:
       outdir=working_dir+"/trimmed"
       
    threads: 32
    
    conda:
        "envs/trim_galore.yml"
       
    shell:
        "mkdir -p {params.outdir} &&"
        "trim_galore --paired {input.R1} {input.R2} --fastqc --output_dir {params.outdir} --cores {threads} "
        ">{log} 2>&1"


        


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
        R1=working_dir+"/trimmed/{sample}_R1_001_val_1.fq.gz",
        R2=working_dir+"/trimmed/{sample}_R2_001_val_2.fq.gz",
        genome=working_dir+"/bwa_index/"+ref_genome_name

    output:
        file=working_dir+"/bwa/{sample}.bam",
        stats=working_dir+"/bwa/{sample}.txt"
        
    log:
        log_dir+"/bwa/{sample}.log"
        
    params:
        outdir=working_dir+"/bwa"
        
    threads: 32
        
    conda:
        "envs/bwa.yml"

    shell:
        "mkdir -p {params.outdir} &&"
        #"bwa-mem2 index {input.genome}&& "
        "bwa-mem2 mem -M -t {threads} {input.genome} {input.R1} {input.R2} > {output.file} "
        "2>{log} && "
        "samtools stats {output.file} >{output.stats}"

   
rule samtools:
    input:
        working_dir+"/bwa/{sample}.bam",

    output:
        file=working_dir+"/samtools/{sample}.bam",
        
    log:
        log_dir+"/samtools/{sample}.log"
        
    params:
        outdir=working_dir+"/samtools"
        
    threads: 32
        
    conda:
        "envs/samtools.yml"  

    shell:
        "mkdir -p {params.outdir} &&"
        "samtools sort --threads {threads} -o {output.file} {input} "
        ">{log} 2>&1"
        
        
rule picard:
    input:
        working_dir+"/samtools/{sample}.bam",

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
        "picard MarkDuplicates REMOVE_DUPLICATES=true I={input} O={output.file} M={output.metrics} CREATE_INDEX=true "
        ">{log} 2>&1"
    
   
rule multiqc:
    input: 
        expand(working_dir+"/fastqc/{sample}_R1_001_fastqc.zip",
               sample=SAMPLES),
        expand(working_dir+"/fastqc/{sample}_R2_001_fastqc.zip",
               sample=SAMPLES),
               
        expand(working_dir+"/trimmed/{sample}_R1_001_val_1_fastqc.zip",
               sample=SAMPLES),
        expand(working_dir+"/trimmed/{sample}_R2_001_val_2_fastqc.zip",
               sample=SAMPLES),
               
        expand(working_dir+"/trimmed/{sample}_R1_001.fastq.gz_trimming_report.txt",
               sample=SAMPLES),
        expand(working_dir+"/trimmed/{sample}_R2_001.fastq.gz_trimming_report.txt",
               sample=SAMPLES),
               
        expand(log_dir+"/trimmed/{sample}.log",
               sample=SAMPLES),
        expand(log_dir+"/trimmed/{sample}.log",
               sample=SAMPLES),
               
        expand(working_dir+"/bwa/{sample}.txt",
               sample=SAMPLES),    
        expand(working_dir+"/picard/{sample}_deduplicated.metrics.txt",
               sample=SAMPLES),

     
    output:
        config['output_dir']+"/multiqc/multiqc.html"
        
    log:
        log_dir+"/multiqc/multiqc.log"
        
    params:
        outdir=config['output_dir']+"/multiqc"

    conda:
        "envs/multiqc.yml"
        
    shell:
        "mkdir -p {params.outdir} &&"
        "multiqc --config multiqc.yaml -o {params.outdir} -n multiqc.html {input} "
        ">{log} 2>&1"

        

