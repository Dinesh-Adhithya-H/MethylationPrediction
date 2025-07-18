configfile: "config/config.yaml"

# Import necessary modules
import os
import subprocess

OUTPUT_DIR=config['OUTPUT_DIR']
FILES_DIR=config['FILES_DIR']

sample_download_links = [i.rstrip('\n').replace('ftp:/','http://') for i in  open(config['TXT_FILE_DIR']).readlines()]
samples = [i.split("/")[-1].rstrip('\n').rstrip(".cram").rstrip(".bam") for i in  open(config['TXT_FILE_DIR']).readlines()]


# Define the rule 'get_prediction' for making methylation state predictions
rule get_model:
    input:
        # Input Python script for prediction
        training_py=config['constants']['PYTHON_DIR']+"MachineLearning_Training.py",
        # Input final CSV file for each sample
        finalfile_dir=OUTPUT_DIR+"combined_final_file.bed",
        # Methylation annotation file
        annotation_file=config['METHYLATION_ANNOTATION_DIR']
    output:
        # Input trained model file
        model_dir=config['MODEL_DIR']
    shell:
        # Execute the Python script to make methylation predictions using the model
        "python {input.training_py} {output.model_dir} {input.finalfile_dir} {input.annotation_file}"

# Define the rule 'combine_final_file' for combining combine_final_file outputs for all samples, to create a pooled dataset.

rule combine_final_file:
    input:
        expand(OUTPUT_DIR+"{sample}/final_file.bed", sample=samples)
    output:
        combined_final_file=OUTPUT_DIR+"combined_final_file.bed"
    run:
        import pandas as pd
        df_combined = pd.DataFrame()
        for i in input:
            data=pd.read_csv(i,low_memory=False,sep="\t")
            df_combined = pd.concat([data,df_combined])
        df_combined = df_combined.groupby(by=["chr","start","end"]).mean()
        df_combined.to_csv(output.combined_final_file,index=True,sep="\t")

# Define the rule 'download_bam_or_cram' for downloading BAM or CRAM files

rule download_bam_or_cram:
    params: link = lambda wildcards: sample_download_links[samples.index(wildcards.sample)],
            fasta_file = config['FASTA_FILE_DIR']
    output: FILES_DIR+"{sample}.bam"
    run:
        if params.link.endswith(".cram"):
            cram_file = FILES_DIR+wildcards.sample+".cram"
            if os.path.exists(params.link):
                shell("samtools view -b -T {params.fasta_file} -o  {output} {cram_file}")
                shell("rm -f {cram_file}")
            else:
                shell("wget {params.link} -O {cram_file}")
                shell("samtools view -b -T {params.fasta_file} -o  {output} {cram_file}")
                shell("rm -f {cram_file}")
        elif params.link.endswith(".bam"):
            if os.path.exists(params.link)==False:
                shell("wget {params.link} -O {output}")
            else:
                pass
        else:
            raise Exception("File format not supported, please give only bam or cram files")


# Define the rule 'bam_process' for processing BAM files
rule bam_process:
    input:
        # Input BAM file for each sample
        bam_file=FILES_DIR+"{sample}.bam",
        # Additional input files required for processing
        CpG_isl_bed_file=config['CpG_ISL_BED_FILE_DIR'],
        fasta_file=config['FASTA_FILE_DIR'],
        remove_reads=OUTPUT_DIR+'{sample}/chromosome_ranges.bed' #config['REMOVE_READS_DIR']

    output: 
        # Output files for each sample after processing
        output_reads=OUTPUT_DIR+"{sample}/output_reads.bed",
        output_reads2=OUTPUT_DIR+"{sample}/output_reads2.bed",
        output_bed=OUTPUT_DIR+"{sample}/output.bed",
        output_fasta=OUTPUT_DIR+"{sample}/output.fasta",
        filtered_output_bed=OUTPUT_DIR+"{sample}/filtered_output.bed",
        output3=OUTPUT_DIR+"{sample}/output3.bed"
    shell:
        """
        # Step 1: Extract reads from BAM file and save as output.bed
        ####### samtools view -f 2 -F 3868 {input.bam_file} | awk '{{print $3 "\t"  $4-2 "\t" $4}}' > {output.output_bed}
        samtools view -F 3868 {input.bam_file} | awk 'BEGIN{OFS="\t"} /^@/ {print; next} $6 !~ /^[0-9]+[SH]/ {print}' | awk '{{print $3 "\t"  $4-2 "\t" $4}}' > {output.output_bed}
        
        # Step 2: Add chr to the read's chromosome
        awk -F'\t' -v OFS='\t' '$1 !~ /^chr/ {{ $1 = "chr" $1 }} 1' {output.output_bed} > {output.filtered_output_bed}

        # Step 3: Remove reads in remove_reads.bed from output.bed and save as output3.bed
        grep -vFf {input.remove_reads} {output.filtered_output_bed} > {output.output3}

        # Step 4: Perform bedtools intersect with CpG island bed file and save as output_reads.bed
        bedtools intersect -b {input.CpG_isl_bed_file} -a {output.output3} -wa -wb > {output.output_reads}

        # Step 5: Perform bedtools intersect with CpG island bed file and save as output_reads2.bed
        bedtools intersect -b {input.CpG_isl_bed_file} -a {output.output3} -wa > {output.output_reads2}

        # Step 6: Extract fasta sequences using bedtools getfasta and save as output.fasta
        bedtools getfasta -fi {input.fasta_file} -bed {output.output_reads2} > {output.output_fasta}
        """

rule extract_chromosomes:
    input:
        bam_file=FILES_DIR+"{sample}.bam"  # Replace with your BAM file
    output:
        OUTPUT_DIR+"{sample}/chromosome_ranges.bed"
    shell:
         """
        samtools view -H {input.bam_file} | grep "^@SQ" | cut -f 2 -d ':' | cut -f 2 -d '@' | cut -f 1 | while read -r chromosome; do
            if [[ "$chromosome" != chr* ]]; then
                chromosome="chr$chromosome"
            fi
            echo -e "$chromosome\\t-1\\t1" >> {output}
        done
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
        python_code=config['constants']['PYTHON_DIR']+"FeatureExtraction.py",
        # Input fasta data directory for each sample
        dinucletide_data_dir=OUTPUT_DIR+"{sample}/output.fasta",
        # Additional input files required for processing
        cpg_island_data_dir=OUTPUT_DIR+'output_cpg.fasta',
        overlap_cpg_dinucleotide_dir=OUTPUT_DIR+"{sample}/output_reads.bed"
    output:
        # Output CSV file for each sample after processing
        OUTPUT_DIR+"{sample}/final_file.bed"
    params:
        # Parameter specifying output directory for each sample
        output_dir_name=OUTPUT_DIR+"{sample}/"
    shell:
        # Execute the Python script to process fasta data and generate final CSV file
        "python {input.python_code} {input.dinucletide_data_dir} {input.cpg_island_data_dir} {input.overlap_cpg_dinucleotide_dir} {params.output_dir_name}"


rule download_fasta_file:
    output: fasta = config['FASTA_FILE_DIR'], index = config['FASTA_FILE_DIR']+".fai"
    params: link = config['FASTA_FILE_LINK'], index_link = config['FASTA_FILE_INDEX_LINK']
    run:
        if os.path.exists(output.fasta)==False:
            shell("wget -O {output.fasta} {params.link}")
        if os.path.exists(output.index)==False:
            if params.index_link!="None":
                shell("wget -O {output.index} {params.index_link}")
            else:
                shell("samtools faidx {output.fasta}")
