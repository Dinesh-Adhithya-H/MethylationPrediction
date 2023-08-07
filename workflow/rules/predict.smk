include: "base.smk"


# Define the rule 'all', which specifies the final output files to be generated
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
        expand(OUTPUT_DIR+"{sample}/methylation_outputs.csv", sample=samples),
        OUTPUT_DIR+"combined_final_file.csv"

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

# Define the rule 'combine_final_file' for combining combine_final_file outputs for all samples, to create a pooled dataset.

rule combine_final_file:
    input:
        expand(OUTPUT_DIR+"{sample}/final_file.csv", sample=samples)
    output:
        combined_final_file=OUTPUT_DIR+"combined_final_file.csv"
    run:
        import pandas as pd
        df_combined = pd.DataFrame(columns=['chr','start','end'])
        for i in input:
            data=pd.read_csv(i,low_memory=False)
            df_combined = pd.merge(data,df_combined,on=['chr','start','end'], how="outer")
        df_combined.to_csv(output.combined_methylation_outputs,index=False)