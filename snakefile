# Snakemake workflow for processing BAM files and predicting methylation state of CpG islands.

singularity: "docker://continuumio/miniconda3:4.7.12"

# Load the configuration settings from the "config.yaml" file

configfile: "config.yaml"

# Import necessary modules
import os

sample_download_links = [i.rstrip('\n').replace('ftp:/','http://') for i in  open(config['TXT_FILE_DIR']).readlines()]
samples = [i.split("/")[-1].rstrip('\n').rstrip(".cram") for i in  open(config['TXT_FILE_DIR']).readlines()]


OUTPUT_DIR=config['OUTPUT_DIR']
FILES_DIR=config['FILES_DIR']

# Define the rule 'all', which specifies the final output files to be generated
rule all:
    input:
        # Expand function generates output file paths for each sample in 'samples'
        #
        expand(FILES_DIR+"{sample}.bam", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output_reads.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output_reads2.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output.fasta", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output3.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/final_file.csv", sample=samples),
        expand(OUTPUT_DIR+"{sample}/methylation_outputs.csv", sample=samples)
        
rule download_cram:
    params:
        link = lambda wildcards: sample_download_links[samples.index(wildcards.sample)]
    output:
        FILES_DIR+"{sample}.cram"
    shell:
        "wget {params.link} -O {output}"

rule convert_cram_to_bam:
    input:
        cram=FILES_DIR+"{sample}.cram",
        fasta_file=config['HOME_DIR']+config['FASTA_FILE_DIR']
    output:
        FILES_DIR+"{sample}.bam"
    params:
        threads=10
    shell:
        "samtools view -bS -@ {params.threads}  -T {input.fasta_file} -o  {output} {input.cram}"

# Define the rule 'bam_process' for processing BAM files
rule bam_process:
    input:
        # Input BAM file for each sample
        bam_file_dir=FILES_DIR+"{sample}.bam",
        # Additional input files required for processing
        CpG_isl_bed_file=config['CpG_ISL_BED_FILE_DIR'],
        fasta_file=config['FASTA_FILE_DIR'],
        remove_reads=config['REMOVE_READS_DIR']
    output: 
        # Output files for each sample after processing
        output_reads=OUTPUT_DIR+"{sample}/output_reads.bed",
        output_reads2=OUTPUT_DIR+"{sample}/output_reads2.bed",
        output_bed=OUTPUT_DIR+"{sample}/output.bed",
        filtered_output_bed=OUTPUT_DIR+"{sample}/filtered_output.bed",
        output_fasta=OUTPUT_DIR+"{sample}/output.fasta",
        output3=OUTPUT_DIR+"{sample}/output3.bed"
    shell:
        """
        # Step 1: Extract reads from BAM file and save as output.bed
        samtools view -f 2 -F 3868 {input.bam_file_dir} | awk '{{print $3 "\t"  $4-2 "\t" $4}}' > {output.output_bed}
        
        # Step 2: Extract reads which belong from chr1 to chr23
        awk '($1 ~ /^chr([1-9]|1[0-9]|2[0-3])$/)' {output.output_bed} > {output.filtered_output_bed}

        # Step 3: Remove reads in remove_reads.bed from output.bed and save as output3.bed
        grep -vFf {input.remove_reads} {output.filtered_output_bed} > {output.output3}

        # Step 4: Perform bedtools intersect with CpG island bed file and save as output_reads.bed
        bedtools intersect -b {input.CpG_isl_bed_file} -a {output.output3} -wa -wb > {output.output_reads}

        # Step 5: Perform bedtools intersect with CpG island bed file and save as output_reads2.bed
        bedtools intersect -b {input.CpG_isl_bed_file} -a {output.output3} -wa > {output.output_reads2}

        # Step 6: Extract fasta sequences using bedtools getfasta and save as output.fasta
        bedtools getfasta -fi {input.fasta_file} -bed {output.output_reads2} > {output.output_fasta}
        """


rule generate_output_cpg:
    input:
        CpG_isl_bed_file=config['CpG_ISL_BED_FILE_DIR'],
        fasta_file=config['FASTA_FILE_DIR'],
    output:
        output_cpg=OUTPUT_DIR+'output_cpg.fasta'
    shell:
        'bedtools getfasta -fi {input.fasta_file}  -bed {input.CpG_isl_bed_file} > {output.output_cpg}'
        

# Define the rule 'bed_process' for further processing of fasta sequences
rule bed_process:
    input:
        # Input Python script for processing
        python_code=config['PYTHON_DIR']+"FeatureExtraction.py",
        # Input fasta data directory for each sample
        dinucletide_data_dir=OUTPUT_DIR+"{sample}/output.fasta",
        # Additional input files required for processing
        cpg_island_data_dir=OUTPUT_DIR+'output_cpg.fasta',
        overlap_cpg_dinucleotide_dir=OUTPUT_DIR+"{sample}/output_reads.bed"
    output:
        # Output CSV file for each sample after processing
        OUTPUT_DIR+"{sample}/final_file.csv"
    params:
        # Parameter specifying output directory for each sample
        output_dir_name=OUTPUT_DIR+"{sample}/"
    shell:
        # Execute the Python script to process fasta data and generate final CSV file
        "python {input.python_code} {input.dinucletide_data_dir} {input.cpg_island_data_dir} {input.overlap_cpg_dinucleotide_dir} {params.output_dir_name}"

# Define the rule 'get_prediction' for making methylation state predictions
rule get_prediction:
    input:
        # Input Python script for prediction
        prediction_py=config['PYTHON_DIR']+"MachineLearning_Prediction.py",
        # Input trained model file
        model_dir=config['MODEL_DIR']+'model.joblib',
        # Input final CSV file for each sample
        finalfile_dir=OUTPUT_DIR+"{sample}/final_file.csv"
    output:
        # Output CSV file with methylation predictions for each sample
        methylation_pred_dir=OUTPUT_DIR+"{sample}/methylation_outputs.csv"
    shell:
        # Execute the Python script to make methylation predictions using the model
        "python {input.prediction_py} {input.model_dir} {input.finalfile_dir} {output.methylation_pred_dir}"