import os
import subprocess
import sys
import shutil
import pyfiglet
import socket
import sendgrid
from termcolor import colored
from concurrent.futures import ThreadPoolExecutor
import time
import gzip
import csv
import pandas as pd
import glob
from pathlib import Path



#-----------------------------------------------------------------------
# ANSI color codes
RED = "\033[31m"  # Red text
GREEN = "\033[32m"
BLUE = "\033[36m"  # Blue text
PURPLE = "\033[35m" # Purple text
RESET = "\033[0m"  # Reset to default color
#-----------------------------------------------------------------------
def run_command(command, retries=3):
    attempt = 0
    while attempt < retries:
        try:
            subprocess.run(command, check=True, shell=True)
            return  # Command succeeded
        except subprocess.CalledProcessError as e:
            attempt += 1
            print(colored(f"An error occurred (attempt {attempt}/{retries}): {e}", "red"), file=sys.stderr)
            if attempt >= retries:
                print(colored(f"Failed to execute command after {retries} attempts. Exiting.", "red"))
                sys.exit(1)
        except KeyboardInterrupt:
            print(colored("\nAnalysis interrupted by user. Exiting.", "red"))
            sys.exit(1)

def check_dependency(command):
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    except subprocess.CalledProcessError:
        print(colored(f"Dependency check failed for command: {command}. Please install the required tool and ensure it is in your PATH.", "red"))
        sys.exit(1)
        
check_dependency("samtools --version")
check_dependency("bowtie2 --version")
check_dependency("fastp --version")


def delete_sbam_files(directory):
    # Get list of files in the directory
    files = os.listdir(directory)
    
    # Iterate through each file
    for file in files:
        # Check if file ends with ".sam"
        if file.endswith(".sam") or file.endswith('.bam') and not file.endswith('.sorted.bam'):
            # Construct the full path of the file
            file_path = os.path.join(directory, file)
            try:
                # Attempt to remove the file
                os.remove(file_path)
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")

def extract_methylation_sites(input_pileup, output_bed):
    with gzip.open(input_pileup, 'rt') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            pileup = fields[4]

            # Check for methylation in the pileup string
            if 'C' in pileup:
                # Assuming a simple BED format with chrom, start, end
                # You might need to adjust this based on your specific requirements
                outfile.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                
def format_bed_file(input_pileup, bed_file_path):
    # Extract methylation sites (assuming this function is already defined)
    extract_methylation_sites(input_pileup, bed_file_path)
    
    # Define the column names for the BED file
    columns = ['chromosome', 'start', 'end', 'name', 'score', 'strand']
    
    # Specify the data types for each column
    dtype_spec = {
        'chromosome': str,
        'start': int,
        'end': int,
        'name': str,
        'score': float,
        'strand': str
    }
    
    # Read the extracted data into a DataFrame with specified dtypes and low_memory=False
    df = pd.read_csv(bed_file_path, sep='\t', header=None, names=columns, dtype=dtype_spec, low_memory=False)
    
    # Remove 'chr' prefix from the chromosome column if present
    df['chromosome'] = df['chromosome'].apply(lambda x: x.replace('chr', '') if isinstance(x, str) and x.startswith('chr') else x)
    
    # Save the formatted DataFrame to the bed file path
    df.to_csv(bed_file_path, sep='\t', header=False, index=False)
    
    print(f"Methylation sites extracted and saved to {bed_file_path}.\n")

def compress_file(input_file, output_file):
    with open(input_file, 'rb') as f_in:
        with gzip.open(output_file, 'wb') as f_out:
            f_out.writelines(f_in)

def bed_to_csv(input_bed, output_csv):
    # Define the awk command to convert BED to CSV and remove 'chr' prefix
    command = f"awk 'BEGIN {{OFS=\",\"}} {{gsub(/^chr/, \"\", $1); print $1, $2, $3, $4, $5, $6}}' {input_bed} > {output_csv}"
    
    # Run the command using run_command function
    result = run_command(command)
    print(result)
     

#----------------------------------------------------------------------

def print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path):
    for chromosome in chromosomes_list:
        bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa" if chromosome != 'hg38' else f"{bwa_base_path}hg38/GRCh38_reference.fa"
        bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie" if chromosome != 'hg38' else f"{bowtie_base_path}hg38/bowtie"
        print(colored(f"\nPaths for chromosome {chromosome}:", "cyan"))
        print(colored("BWA Chromosome Path:", "cyan"), bwa_chrom_path)
        print(colored("Bowtie Index Path:", "cyan"), bowtie_index_path)

def read_accession_numbers(file_path):
    try:
        with open(file_path, 'r') as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        print(colored("The specified file was not found. Please check the file path and try again.", "red"))
        sys.exit(1)
    except Exception as e:
        print(colored(f"An error occurred while reading the file: {e}", "red"))
        sys.exit(1)

def is_file_empty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) <= 800

def delete_intermediate_files(accession_number, chromosome):
    intermediate_files = [
        f"{accession_number}/{accession_number}.fastq",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.sam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf",
        f"{accession_number}/{accession_number}_fastqc.zip"
        f"{accession_number}/{accession_number}_mapped_{chromosome}.bam",
    ]
    for file_path in intermediate_files:
        if os.path.isfile(file_path):
            os.remove(file_path)
            print(colored(f"Deleted {file_path}", "green"))

def clear_directory(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)

def get_verified_path(prompt_message):
    while True:
        path = input(colored(prompt_message, "cyan")).strip()
        if os.path.exists(path):
            return path
        else:
            print(colored("The provided path does not exist. Please try again.", "yellow"))

def sort_bam_file(bam_file_path, output_file_path, threads=4):
    command = [
        "samtools", "sort",
        "-o", output_file_path,
        "-@", str(threads),  # Specify the number of threads
        bam_file_path
    ]
    run_command(" ".join(command))

def parallel_sort_bam_files(bam_files, threads_per_sort=2):
    """Sort multiple BAM files in parallel."""
    with ThreadPoolExecutor() as executor:
        futures = []
        for bam_file in bam_files:
            output_file_path = bam_file.replace(".bam", ".sorted.bam")
            # Submit a sorting task for each BAM file
            futures.append(executor.submit(sort_bam_file, bam_file, output_file_path, threads_per_sort))
        
        # Wait for all futures to complete (optional, depending on your use case)
        for future in futures:
            future.result()  # This will re-raise any exception that occurred in the thread

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

text = "V&S METHYL"
font = "banner3-D"
terminal_width = os.get_terminal_size().columns
f = pyfiglet.Figlet(font=font, width=terminal_width)
logo = f.renderText(text)
print(logo.center(terminal_width))
final_bam_file = ''
#methylation_file_path = '/media/cancer10/VictorUSB/Methylation/methylation_results3/result4/methylation_extration_v1.5.py' :)-Victor Wang
bwa_base_path = "/usr/local/bin/bwa/"
bowtie_base_path = "/usr/local/bin/bowtie/"
fastp_path = "/usr/bin/fastp"

if not os.path.exists(bwa_base_path):
    bwa_base_path = get_verified_path("1. BWA base path not found. Please enter the correct BWA base path: ")
if not os.path.exists(bowtie_base_path):
    bowtie_base_path = get_verified_path("2. Bowtie base path not found. Please enter the correct Bowtie base path: ")


print(colored("1. Please enter the path to the accession list file: ", "white"))
accession_list_file = input().strip()

print(colored("2. Are the reads single-end (1) or paired-end (2)? Enter 1 or 2: ", "white"))
read_type = input().strip()

accession_numbers = read_accession_numbers(accession_list_file)
total_accession_numbers = len(accession_numbers)
all_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

print(colored(f"\nTotal accession numbers in the file: {total_accession_numbers}", "cyan"))

chromosomes_input = input(colored("3. Please enter the chromosomes to be analyzed, separated by a comma, or type 'all' to analyze all chromosomes: ", "white"))
chromosomes_list = all_chromosomes if chromosomes_input.lower() == 'all' else [chromosome.strip() for chromosome in chromosomes_input.split(',')]

print(colored("4. Are the SRR files normal or canerous? n/c ", "white"))
methylation_naming = input().strip()
normal = False
if methylation_naming == "n":
    normal == True

print(colored("5. Do you want to run the methlyation program? y/n ", "white"))
methylation_running = input().strip()

accession_numbers_to_analyze = []
start_time = time.time()
current_directory = os.getcwd()

for accession_number in accession_numbers:
    print_information = True
    methylation_path = os.path.join(current_directory, accession_number, "hello.bed") #remember to change back
    #methylation_path = os.path.join(current_directory, accession_number, accession_number + ".bed") #remember to change back

    if os.path.isfile(methylation_path) and not is_file_empty(methylation_path):
        print(colored(f"Methylation file for {accession_number} already exists. Skipping analysis for this accession number...", "yellow"))  
    else:
        accession_numbers_to_analyze.append(accession_number)
        if print_information == True:
            print(colored(f"Final Methylation file in {accession_number} is missing.", "green"))
            print_information = False

print(colored(f"Accession numbers pending analysis: {len(accession_numbers_to_analyze)}", "green"))

if not accession_numbers_to_analyze:
    print(colored("No accession numbers left to analyze.", "red"))

    print(colored("Do you want to run it again? y/n", "white"))
    continue_input = input().strip()
    if continue_input == 'y':
        for accession_number in accession_numbers:
            accession_numbers_to_analyze.append(accession_number)
    else:
        sys.exit(0)

num_to_analyze = int(input(colored("6. How many accession numbers do you want to analyze? ", "white")))
accession_numbers_to_analyze = accession_numbers_to_analyze[:num_to_analyze]

clear_directories = [d for d in os.listdir(current_directory) 
                     if os.path.isdir(os.path.join(current_directory, d)) 
                     and (d.startswith("SRR") or d.startswith("ERR"))]

print("Directories to clear:", clear_directories)

if clear_directories:
    for c in clear_directories:
        c_path = os.path.join(current_directory, c)
        # Correctly list files in the subdirectory
        files_need_cleared = [f for f in os.listdir(c_path) if "tmp" in f or "sam" in f]
        
        for f in files_need_cleared:
            clear_path = os.path.join(c_path, f)
            try:
                # Attempt to remove the file
                os.remove(clear_path)
                print(f"{clear_path} cleared :)")
            except Exception as e:
                print(f"Error deleting {clear_path}: {e}")
        

print(colored("List of chromosomes to be analyzed:", "cyan"), chromosomes_list)
print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path)
for accession_number in accession_numbers_to_analyze:
    # Define file paths for single-end and paired-end reads
    trimmed_file_single = f"{accession_number}/{accession_number}_trimmed.fq.gz"
    trimmed_file_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq.gz"
    trimmed_file_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq.gz"
    trimmed_file_path = os.path.join(current_directory,accession_number,trimmed_file_paired_2)
    if not os.path.isfile(trimmed_file_path):
        # Download and trimming logic based on read type
        if read_type == '1' and not os.path.isfile(trimmed_file_single):
            print(colored(f"\nDownloading single-end sequence number {accession_number} from SRA...", "cyan"))
            run_command(f"fastq-dump {accession_number} -O {accession_number}/")
            print(colored(f"\nTrimming {accession_number} with fastp...", "cyan"))
            trim_command = f"{fastp_path} -i {accession_number}/{accession_number}.fastq -o {trimmed_file_single} --thread=4"
            run_command(trim_command)
            shutil.move("fastp.html", f"{accession_number}/fastp.html")
            shutil.move("fastp.json", f"{accession_number}/fastp.json")
        elif read_type == '2' and not (os.path.isfile(trimmed_file_paired_1) and os.path.isfile(trimmed_file_paired_2)):
            print(colored(f"\nDownloading paired-end sequence number {accession_number} from SRA...", "cyan"))
            run_command(f"fastq-dump --split-files {accession_number} -O {accession_number}/")
            print(colored(f"Trimming {accession_number} with fastp...", "cyan"))
            trim_command = f"{fastp_path} -i {accession_number}/{accession_number}_1.fastq -I {accession_number}/{accession_number}_2.fastq -o {trimmed_file_paired_1} -O {trimmed_file_paired_2} --thread=4"
            run_command(trim_command)
            shutil.move("fastp.html", f"{accession_number}/fastp.html")
            shutil.move("fastp.json", f"{accession_number}/fastp.json")
        else:
            print(colored(f"Trimmed files for {accession_number} found or invalid read type. Skipping download and trimming.", "yellow"))

        for chromosome in chromosomes_list:
            final_bam_file = f"{accession_number}/{accession_number}_mapped_{chromosome}.sorted.bam"
            if os.path.isfile(final_bam_file) and not is_file_empty(final_bam_file):
                print(colored(f"Bam file for {accession_number}, chromosome {chromosome} already exists. Skipping analysis...", "yellow"))
                continue
            elif is_file_empty(final_bam_file):
                print(colored(f"Bam file for {accession_number}, chromosome {chromosome} is empty. Deleting and adding to analysis...", "yellow"))
                os.remove(final_bam_file)
            else:
                # Set paths for BWA and Bowtie2 indices
                bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa" if chromosome != 'hg38' else f"{bwa_base_path}hg38/GRCh38_reference.fa"
                bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie" if chromosome != 'hg38' else f"{bowtie_base_path}hg38/bowtie"

                # Set the correct input file for mapping based on read type
                input_file_for_mapping = trimmed_file_single if read_type == '1' else trimmed_file_paired_1

                print(colored(f"\nMapping {accession_number} reads using Bowtie2 for chromosome {chromosome}...", "cyan"))

                bowtie_command = f"bowtie2 --very-fast-local -x {bowtie_index_path} -U {input_file_for_mapping} -S {accession_number}/{accession_number}_mapped_{chromosome}.sam"
                run_command(bowtie_command)

                run_command(f"samtools view -S -b {accession_number}/{accession_number}_mapped_{chromosome}.sam > {accession_number}/{accession_number}_mapped_{chromosome}.bam")

                print(colored("\nSorting using Samtools...", "magenta"))
        
        for dir_name in os.listdir(current_directory):
            dir_path = os.path.join(current_directory, dir_name)
            if os.path.isdir(dir_path) and dir_name.startswith("SRR"):
                print(f"Processing directory: {dir_name}")
                os.chdir(dir_path)  # Change working directory to the SRR directory

                # Initialize a list to hold .bam files that need sorting
                bam_files_to_sort = []

                # Check each .bam file to see if a sorted version already exists
                for f in os.listdir(dir_path):
                    if f.endswith(".bam") and not f.endswith(".sorted.bam"):
                        sorted_version = f.replace(".bam", ".sorted.bam")
                        if not os.path.exists(sorted_version):  # Check if sorted file does not exist
                            bam_files_to_sort.append(f)

                # Run parallel sorting if there are BAM files to sort
                if bam_files_to_sort:
                    parallel_sort_bam_files(bam_files_to_sort, threads_per_sort=2)
                else:
                    print("No unsorted .bam files to sort in this directory.")

                os.chdir(current_directory)  # Change back to the original directory

            # Delete intermediate files to save disk space
            delete_intermediate_files(accession_number, chromosome)

if methylation_running == "n":
    end_time = time.time()
    elapsed_time_seconds = end_time - start_time
    hours = int(elapsed_time_seconds // 3600)
    minutes = int((elapsed_time_seconds % 3600) // 60)
    print(colored(f'Program finished in {hours}hours and {minutes}minutes','red'))

if methylation_running == "y":
    print("Tring to run the methylation program...")
    current_directory = os.getcwd()

    # Get a list of directories starting with SRR or ERR in the current directory
    directories = [d for d in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, d)) and (d.startswith("SRR") or d.startswith("ERR"))]
    # Iterate over the directories
    for directory in accession_numbers_to_analyze:
        # cd into the working directory
        os.chdir(str(directory))
        print(GREEN + f'Currently working with {str(directory)}\n' + RESET)
        bam_directory = os.path.join(current_directory,directory,"bam")
        print(BLUE + f'Processing the bismark_bowtie2 program in {str(directory)}... \n' + RESET)
        input_directory = os.path.join(current_directory, directory)
        
        
        ##### Bismark_Bowtie2 #####
        
        # Check whether bam directory has already been generated (avoid repetitive operation)
        if not os.path.isdir(bam_directory):
            # Build the full path of the input directory
            os.makedirs(bam_directory, exist_ok=True)
            # Get a list of files ending with _trimmed.fq.gz in each directory
            files = [f for f in os.listdir(input_directory) if f.endswith("_trimmed.fq.gz")]
    
            # Iterate over the files in the directory
            for file in files:
                    # Build the command to run Bismark in the specific directory
                    command_bowtie2 = [
                        "bismark",
                        "--bowtie2",
                        "/usr/local/bin/bowtie/hg38_bismark/",
                        "-o",
                        "bam",
                        os.path.join(input_directory, file)
            ]
                    # Run the bowtie2
                    run_command(command_bowtie2)
                    print('\n')
        print(BLUE + f'Bismark_bowtie2 program in {str(directory)} is done or has already been processed \n' + RESET)
        
        
        ##### Samtool_sorting #####
        '''
        print(BLUE + f'Processing the samtool_sorting program in {str(directory)}... \n' + RESET)
        # Check whether 'sorted_alignment.bam' has already been generated (avoid repetitive operation)
        output_file_path = os.path.join(bam_directory,'sorted_alignment.bam')
        
        if os.path.exists(bam_directory) and os.path.isdir(bam_directory) and not os.path.exists(output_file_path):
            # Get a list of files ending with .bam in the "bam" directory
            bam_files = [f for f in os.listdir(bam_directory) if f.endswith(".bam")]

            # Iterate over the files in the "bam" directory
            for bam_file in bam_files:
                # Build the command to run samtools sort
                command_samtool_sort = [
                    "samtools",
                    "sort",
                    "-o",
                    os.path.join(bam_directory, "sorted_alignment.bam"),
                    os.path.join(bam_directory, bam_file)
                ]

                # Run the samtool
                subprocess.run(command_samtool_sort)
                print('\n')
        print(BLUE + f'Samtool_sorting program in {str(directory)} is done or has already been processed \n' + RESET)   
        '''
        ##### Samtool_merging #####
        
         # Iterate through all items in the current_director
        final_merged_file = os.path.join(bam_directory,"sorted_alignment.bam")
        if not os.path.isfile(final_merged_file):
                # Look for sorted BAM files to merge
            bam_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith(".sorted.bam")]    
            # Ensure there are files to merge
            if bam_files:
                # Merge the files using samtools
                merge_command = ['samtools', 'merge', final_merged_file] + bam_files
                run_command(" ".join(merge_command))  # Convert list to string
                print(f"Merged file created at {final_merged_file}")
            else:
                print("No .sorted.bam files found for merging.")
        else:
            print(f"Final merged file already exists at {final_merged_file}. Skipping merging.")


        ##### Samtool_pileup #####
        bam_files_sorted = [f for f in os.listdir(bam_directory) if f.endswith("sorted_alignment.bam")]
        print(BLUE + f"Running Samtool_pileup program in {directory}\n" + RESET)
        # Iterate over the sorted_alignment.bam files in the "bam" directory
        pileup_file_path = os.path.join(current_directory,directory,"output.pileup")
        # If there is output.pileup, skip this directory
        if not os.path.exists(pileup_file_path) or is_file_empty(pileup_file_path):
            if is_file_empty(pileup_file_path):
                os.remove(pileup_file_path)
            for bam_file_sorted in bam_files_sorted:
                # Build the command to run samtools mpileup
                command = [
                "samtools",
                "mpileup",
                "-f",
                "/media/cancer10/VictorUSB/hg38/GRCh38_reference.fa",
                os.path.join(bam_directory, bam_file_sorted),
                    ">",
                os.path.join(input_directory, "output.pileup")
                ]

                # Run the command
                run_command(" ".join(command))
        print(BLUE + f"mpileup completed!.\n" + RESET)
                
        ##### Gizp_compression #####
        
        print(BLUE + f"Compressing pileup files... \n" + RESET)
        gzip_file = "output.pileup.gz"
        input_file = "output.pileup"
        #if there is compressed gzip file, skip
        if not os.path.exists(gzip_file):
            #run gzip compression command
            compress_file(input_file, gzip_file)
        print(BLUE + f"Compressing complete \n" + RESET)
        
        
        ##### Methylation_location_extractor #####

        print(BLUE + f"Extracting methylation sites...\n" + RESET)
        
        file_to_delete = os.path.join(current_directory,directory, 'methylation_site.bed')
        if os.path.exists(file_to_delete):
            # Delete the file
            os.remove(file_to_delete)
            print(f'{file_to_delete} has been deleted successfully.')
    
        bed_file_path = os.path.join(current_directory,directory, directory + ".bed")
        #if there is bed file, skip
        if not os.path.exists(bed_file_path):
            input_pileup = os.path.join(current_directory,directory,"output.pileup.gz")
            #run bed file conversion command
            print("Converting .pileup file to .bed file...")
            
            format_bed_file(input_pileup, bed_file_path)

            print(f"Methylation sites extracted and saved to {bed_file_path}.\n")
        
        output_file_path = os.path.join(current_directory,directory, directory + "_chr.bed")
        if not os.path.exists(output_file_path):
            with open(bed_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
                for line in infile:
                    columns = line.strip().split('\t')
                    columns[0] = 'chr' + columns[0]
                    outfile.write('\t'.join(columns) + '\n')
            print(output_file_path + 'has been created for IGV')
        
        ##### Converting .bed to .csv #####
        
        print(BLUE + f"Generating .csv file... \n" + RESET)
        
        csv_file  = [d for d in os.listdir(input_directory) if d.endswith("output.csv")]
        if len(csv_file) == 0:
            
            if normal == True:
                output_csv_file = 'normal_' + str(directory) + '_output.csv'
            else:
                output_csv_file = 'cancer_' + str(directory) + '_output.csv'
            

            bed_to_csv(bed_file_path, output_csv_file)
            print(BLUE + f"Bed file converted to Csv file and saved to {output_csv_file}" + RESET)
            
            file_path = Path(current_directory,directory,output_csv_file)
            
            target_path = file_path.parent.parent / file_path.name
            shutil.move(str(file_path), str(target_path))

        delete_sbam_files(input_directory)                
        os.chdir(str(current_directory))
        print( GREEN + "Now changing into another directory\n" + RESET)
        
    end_time = time.time()
    elapsed_time_seconds = end_time - start_time
    hours = int(elapsed_time_seconds // 3600)
    minutes = int((elapsed_time_seconds % 3600) // 60)
    print(colored(f'Program finished in {hours}hours and {minutes}minutes','red'))

# Define the root directory where directories starting with "SRR" are located
root_dir = Path.cwd()  # Adjust this path as necessary
# Define the path for the directory where all the Parquet files will be stored
parquet_dir = os.path.join(root_dir, 'all_parquet')

# Create the 'all_parquet' directory if it doesn't exist
if not os.path.exists(parquet_dir):
    os.makedirs(parquet_dir)

for filename in accession_numbers_to_analyze:
    if filename.endswith(".csv"):
        # Construct full path of the CSV file
        csv_path = os.path.join(root_dir, filename)
        # Read the CSV file
        try:
            df = pd.read_csv(csv_path)
            # Define the new path with the .parquet extension
            parquet_filename = filename.replace(".csv", ".parquet")
            parquet_path = os.path.join(parquet_dir, parquet_filename)
            # Write the DataFrame to a Parquet file
            df.to_parquet(parquet_path)
            print(f"Converted {csv_path} to {parquet_path}")
        except Exception as e:
            print(f"Failed to convert {csv_path} due to {e}")

print("Conversion completed.")
print( PURPLE + "---ALL processing completed :) ---\n" + RESET)
            
