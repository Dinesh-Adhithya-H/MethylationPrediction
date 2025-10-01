# WGS2meth

## Getting Started

WGS2meth predicts the methylation state of CpG islands directly from standard DNA-seq alignemnt data. It exploits fragmentation biases from library preparation step: methylated CpG dinucleotides are more prone to hydrolisis, which leave characteristic patterns in read start coordinates within CpG islands. By analyzing these patterns, WGS2meth infers whether each CpG island is methylated or unmethylated.

Input: aligned BAM/CRAM files.

## What does it do?
1. WGS2meth takes a .txt file, containing url links or local paths to BAM/CRAM files as an input.
2. If needed, alignemnt files are downloaded and converted to BAM. If no local reference genome is provided, it can download one from a user-supplied url link.
3. Reads are filtered and read positions are extracted from a BAM file.
4. For each of 16 dinucleotides fragmentation rates are measured for each CpG island. Then CpG islands methylation state is predicted using a pre-trained model.
5. Depending on the input settings - aggregation of input files can be made.

This scheme illustrates corresponding snakemake pipeline:

![Workflow](./Figures/dag.svg)

## Installation

1. Clone the GitHub repository to your local machine.
```sh
  git clone https://github.com/Dinesh-Adhithya-H/MethylationPrediction.git
```
2. Go to the directory containing the git clone
```sh
  cd MethylationPrediction
```
3. Run the following command to check if necessary tools such as samtools, bedtools, snakemake, and python are installed.

```sh
bash workflow/envs/check_tools.sh
```
Please, run the bash snippets below to create a virtual Python environment with the necessary modules for the tool.

4. Create a virtual environment
```sh
python -m venv WGS2meth_venv
```
5. Activate the virtual environment
```sh
source WGS2meth_venv/bin/activate
```
6. Install dependencies from python_requirements.txt
```sh
pip install -r workflow/envs/python_requirements.txt
```
7. Check if everything works (snakemake dry run)
```sh
snakemake -n
```

## Usage

1. Edit the config file 'config.yaml', to set up the mode and ensure the right directories and inputs are used.
Some comments on how to set up config file:
``` yaml
CpG_ISL_BED_FILE: 'Bed file with CpG islands coordinates has to be provided. Ones for hg19 and hg38 can be found in the "resources/" directory.'
TXT_FILE_LIST: 'The file which contains paths to the input bam files or urls of bam/cram files to download.'
MODE: 'Could be either "predict" or "train". If "train", the model will be trained and saved in "MODEL". If "predict", the model will be loaded and used to predict the methylation state of input CpG islands.'
METHYLATION_ANNOTATION: 'The csv file which contains the methylation state of the CpG islands. It is needed only in the "train" mode'
RATIO: '"None" or some float value "r" which specifies an expected fraction of methylated CpG islands. It is only used in "predict" mode. If it is set to "None" - the ratio observed during the input model's training will be used.'
MODEL: 'Model from "workflow/models/" folder has to be specified. Depending on the mode, a model is either created ("train") or used as an input ("predict"). Better always keep "model_degault" model intact and overwrite only "model.joblib" file during training.'
SAMPLE_TYPES: 'Whether the samples belong to the same or different sources ("SAME" or "DIFFERENT"). If "SAME" - all outputs are aggregated at the last step in one combined file. Methylation state predctions are also made on the aggregated sample.'
```

2. Edit the 'file_links_local.txt' file with the list of samples (or url links) to run WGS2meth on.

3. Run the snakemake pipeline.
``` sh
  snakemake --cores 10
```

## Outputs

Without aggregation, in "predict" mode there will be two output files:
1. "methylation_outputs.bed", which includes methylation states and ratio "r" adjusted methylation states.
2. "final_file.bed" with dinucleotide's fragmentation rates and frequencies.

## Citation


