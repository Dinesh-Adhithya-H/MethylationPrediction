# Snakefile
configfile: "config/config.yaml"

OUTPUT_DIR=config['OUTPUT_DIR']
FILES_DIR=config['FILES_DIR']

mode = config['MODE']

# Check if the chosen mode is valid
valid_modes = ["train", "predict"]
if mode not in valid_modes:
    raise ValueError(f"Invalid mode selected. Please choose one of: {', '.join(valid_modes)}")

if mode == "train":
    include: "workflow/rules/train.smk"

    OUTPUT_DIR=config['OUTPUT_DIR']
    FILES_DIR=config['FILES_DIR']
    sample_download_links = [i.rstrip('\n').replace('ftp:/','http://') for i in  open(config['TXT_FILE_LIST']).readlines()]
    samples = [i.split("/")[-1].rstrip('\n').rstrip(".cram").rstrip(".bam") for i in  open(config['TXT_FILE_LIST']).readlines()]

    rule all:
        input:
            config['FASTA_FILE'],
            expand(FILES_DIR + "{sample}.bam", sample=samples),
            expand(OUTPUT_DIR + "{sample}/final_file.bed", sample=samples),
            config['MODEL']
        params:
            unnecessary_files = [
            OUTPUT_DIR + 'output_cpg.fasta',
            expand(OUTPUT_DIR + "{sample}/filtered_output.bed", sample=samples),
            expand(OUTPUT_DIR + "{sample}/output_reads.bed", sample=samples),
            expand(OUTPUT_DIR + "{sample}/output_reads2.bed", sample=samples),
            expand(OUTPUT_DIR + "{sample}/output.bed", sample=samples),
            expand(OUTPUT_DIR + "{sample}/output.fasta", sample=samples),
            expand(OUTPUT_DIR + "{sample}/output3.bed", sample=samples),
            expand(OUTPUT_DIR + "{sample}/chromosome_ranges.bed",sample=samples)
            ]
        run:
            shell(f"rm -f {params.unnecessary_files[0]}")
            for files in params.unnecessary_files[1:]:
                for sample_files in files:
                    shell(f"rm -f {sample_files}")


elif mode == "predict":
    include: "workflow/rules/predict.smk"
    
    OUTPUT_DIR=config['OUTPUT_DIR']
    FILES_DIR=config['FILES_DIR']

    sample_download_links = [i.rstrip('\n').replace('ftp:/','http://') for i in  open(config['TXT_FILE_LIST']).readlines()]
    samples = [i.split("/")[-1].rstrip('\n').rstrip(".cram").rstrip(".bam") for i in  open(config['TXT_FILE_LIST']).readlines()]

    if config['SAMPLE_TYPES']=="SAME":
        rule all:
            input:
                config['FASTA_FILE'],
                expand(FILES_DIR+"{sample}.bam", sample=samples),
                expand(OUTPUT_DIR+"{sample}/final_file.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/methylation_outputs.bed", sample=samples),
                OUTPUT_DIR+"combined_final_file.bed",
                OUTPUT_DIR+"combined_methylation_outputs.bed",
                OUTPUT_DIR+"summary_with_combined_final_file.bed"

            params:
                unnecessary_files = [
                OUTPUT_DIR+'output_cpg.fasta',
                expand(OUTPUT_DIR+"{sample}/filtered_output.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output_reads.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output_reads2.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output.fasta", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output3.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/chromosome_ranges.bed",sample=samples)
                ]
            run:
                shell(f"rm -f {params.unnecessary_files[0]}")
                for files in params.unnecessary_files[1:]:
                    for sample_files in files:
                        shell(f"rm -f {sample_files}")


    elif config['SAMPLE_TYPES']=="DIFFERENT":
        rule all:
            input:
                config['FASTA_FILE'],
                expand(FILES_DIR+"{sample}.bam", sample=samples),
                expand(OUTPUT_DIR+"{sample}/final_file.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/methylation_outputs.bed", sample=samples),
                OUTPUT_DIR+"summary.bed"

            params:
                unnecessary_files = [
                OUTPUT_DIR+'output_cpg.fasta',
                expand(OUTPUT_DIR+"{sample}/filtered_output.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output_reads.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output_reads2.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output.fasta", sample=samples),
                expand(OUTPUT_DIR+"{sample}/output3.bed", sample=samples),
                expand(OUTPUT_DIR+"{sample}/chromosome_ranges.bed",sample=samples)
                ]


            run:
                shell(f"rm -f {params.unnecessary_files[0]}")
                for files in params.unnecessary_files[1:]:
                    for sample_files in files:
                        shell(f"rm -f {sample_files}")

    else:
        raise ValueError(f"Invalid sample type selected. Please choose one of: SAME, DIFFERENT")

        
