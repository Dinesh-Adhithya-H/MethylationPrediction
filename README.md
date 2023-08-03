# MethylationPrediction

## Getting Started

Using the information of dinucleotide around which reads break during sequencing, this project aims to predict the methylation state of CpG Islands.

## What does it do?
1. The package takes a .txt file containing links to cram files as input from the 1000 Genomes project.
2. The .cram files are downloaded.
3. Converts the .cram file to a .bam file using a fasta file. The directory of the fasta file shall be edited in the config file.
4. The read positions are extracted from the bam file using bedtools.
5. The nucleotide starting positions are extracted, and then the ratio of the occurring frequency of the dinucleotide to the expected frequency of the dinucleotide.
6. Then a finalfile.csv is generated, with contains the coordinates of the CpG islands and 16 feature values, each corresponding to a dinucleotide.
7. Then a methylation_prediction.csv is generated, where the CpG Island methylation state is predicted using a pre-trained model.

## Installation

''' sh
  git clone https://github.com/Dinesh-Adhithya-H/MethylationPrediction.git
'''

