include: "base.smk"

rule all:
    input:
        # Expand function generates output file paths for each sample in 'samples'
        OUTPUT_DIR+"dependencies_check.txt",
        config['FASTA_FILE_DIR'],
        expand(FILES_DIR+"{sample}.bam", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output_reads.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output_reads2.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output.fasta", sample=samples),
        expand(OUTPUT_DIR+"{sample}/output3.bed", sample=samples),
        expand(OUTPUT_DIR+"{sample}/final_file.csv", sample=samples),
        OUTPUT_DIR+"combined_final_file.csv",
        config["MODEL_FILE"]