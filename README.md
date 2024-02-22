# Methylation_project overview

## Aim
This project is designed to advance our understanding of methylation patterns across genomes, specifically focusing on distinguishing between normal and cancerous methylation. The primary goal is to extract methylation locations from genetic data and utilize a machine-learning model to differentiate between normal and cancerous methylation states.

## Components
The project consists of two main programs and some small programs:
1. **Methylation Location Extraction (unt_trimmed_bin_analysis_bam_v1.9.py)**: A program developed to accurately identify and extract methylation sites from genetic sequences. 
2. **Machine Learning Model for Methylation Analysis (machine_methyl.py)**: A  machine learning algorithm trained to differentiate between normal and cancerous methylation patterns. By learning from a dataset of known methylation states, the model can classify new, unseen methylation sites.

## Language, input, and output of these programs
- We use Python3 as the programming language.
  1. **For unt_trimmed_bin_analysis_bam program**:
- You have to download the accesstion list from NCBI (https://www.ncbi.nlm.nih.gov/sra) that contains the strings starting with "SRR" or "ERR".
- You will get the .csv output, which is a sheet of the location of the methylation
  2. **For machine learning program**:
- You have to transfer .csv files to .parquet files in order to make this program work.
- The output of this learning program is still under construction...
  
## Dependency Installation
**Linux is recommended for these programs to run**
1. **Install anaconda (https://www.anaconda.com/download)**:
2. **Run the command: bash ~{the path to the things you downloaded}**
3. **Run the program dependency_installation.py**
   

## Acknowledgments
- We would like to express our gratitude to Dr. Robert Aguilar for his foundational work on the SNP analysis program, our project builds upon his original code.
- Dr. Robert Aguilar's Github Profile: https://github.com/crisprmax
- Reference code: https://github.com/crisprmax/SNP-identifier-Python

