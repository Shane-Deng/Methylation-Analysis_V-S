from http.client import ImproperConnectionState
import os
import sys
import re
import requests
import sendgrid
import pandas as pd
import pyarrow.parquet as pq
import shutil
from tqdm import tqdm
from termcolor import colored
import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, accuracy_score
from sklearn.utils import shuffle
from imblearn.over_sampling import SMOTE
from sendgrid.helpers.mail import Mail, Email, To, Content

# Function to display the first few rows of a Parquet file
def display_first_rows(file_name, num_rows=5):
    try:
        table = pq.read_table(file_name)
        for row in table.to_pandas().head(num_rows).itertuples(index=False):
            print(row)
        print('\n')
    except Exception as e:
        print(colored(f"Error reading {file_name}: {str(e)}", "red"))

###########################################################################################

# Part 1: Process .vcf files and save as .parquet files

# Scan for all 'comb.vcf' files
vcf_files = [f for f in os.listdir('.') if f.endswith('comb.vcf')]

# Process each file
for file in tqdm(vcf_files, desc="Processing Files"):  # tqdm wraps around the iterable providing a progress bar
    # Read the data
    df = pd.read_csv(file, sep='\t', header=None)
    # Rename columns for clarity
    df.columns = ['Chromosome', 'Position', 'Reference', 'Alternate'] + [f'Individual_{i}' for i in range(1, len(df.columns)-3)]
    # Save or manipulate dataframe as needed
    # Example: Save preprocessed file
    df.to_parquet(f'{file}_procs.parquet')


#############################################################################################

# Part 2: Display the first few rows of .parquet files with the specified pattern
pattern = 'ch*.procs.parquet'
matching_files = [file for file in os.listdir() if file.startswith('ch') and file.endswith('procs.parquet')]

for parquet_file in matching_files:
    print(colored(f"First few rows of {parquet_file}:", "magenta"))
    try:
        table = pq.read_table(parquet_file)
        n = 5
        for row in table.to_pandas().head(n).itertuples(index=False):
            print(row)
        print('\n')
    except Exception as e:
        print(colored(f"Error reading {parquet_file}: {str(e)}", "red"))

#################################################################################################

# Part 3: Automatically detect cancer and normal .parquet files
cancer_files = [f for f in os.listdir() if '_comb.vcf_procs.parquet' in f and 'normal' not in f]
normal_files = [f for f in os.listdir() if 'normal_comb.vcf_procs.parquet' in f]

##################################################################################################

# Part 4: Process each cancer file with the corresponding normal file

# Before processing cancer files, check if they exist
if not cancer_files or not normal_files:
    print(colored("Cancer or normal files missing. Exiting Part 4.", "red"))
    exit()

for cancer_file in tqdm(cancer_files, desc=colored("Filtering Cancer Files", "blue")):
    chromosome, cancer_type = cancer_file.split('_')[:2]
    normal_file = f'{chromosome}_normal_comb.vcf_procs.parquet'

    if normal_file in normal_files:
        # Load data
        cancer_data = pd.read_parquet(cancer_file)
        normal_data = pd.read_parquet(normal_file)

        # Assign labels
        cancer_data['label'] = 1  # Label for cancer data
        normal_data['label'] = 0  # Label for normal data


        # Calculate the frequency of each variant for both cancer and normal data
        cancer_data['Variant_Frequency'] = cancer_data['Individual_1'].apply(lambda x: x.count('1') / len(x))
        normal_data['Variant_Frequency'] = normal_data['Individual_1'].apply(lambda x: x.count('1') / len(x))

        # Create a composite variant identifier for both datasets
        cancer_data['Variant_Key'] = cancer_data['Chromosome'].astype(str) + '_' + cancer_data['Position'].astype(str) + '_' + cancer_data['Reference'] + '_' + cancer_data['Alternate']
        normal_data['Variant_Key'] = normal_data['Chromosome'].astype(str) + '_' + normal_data['Position'].astype(str) + '_' + normal_data['Reference'] + '_' + normal_data['Alternate']

        # Define thresholds for filtering
        common_threshold = 0.05  # Example threshold for too common variants
        rare_threshold = 0.01    # Example threshold for too rare variants

        # Filter out variants too common or too rare in normal cohort
        filtered_normal_data = normal_data[(normal_data['Variant_Frequency'] <= common_threshold) &
                                           (normal_data['Variant_Frequency'] >= rare_threshold)]

        # Filter cancer data based on the filtered variant keys from the normal cohort
        filtered_variants = set(filtered_normal_data['Variant_Key'])
        filtered_cancer_data = cancer_data[cancer_data['Variant_Key'].isin(filtered_variants)]

        # Save the filtered data
        cancer_output_file = f'filtered_{chromosome}_{cancer_type}_cancer.parquet'
        normal_output_file = f'filtered_{chromosome}_normal.parquet'
        filtered_cancer_data.to_parquet(cancer_output_file)
        filtered_normal_data.to_parquet(normal_output_file)


        # Display the first few lines of each filtered file
        print(colored(f"First few rows of {cancer_output_file}:", "green"))
        display_first_rows(cancer_output_file)
        print(colored(f"First few rows of {normal_output_file}:", "green"))
        display_first_rows(normal_output_file)

if os.path.exists(cancer_output_file) and os.path.exists(normal_output_file):
    print(f"Cancer file created: {cancer_output_file}")
    print(f"Normal file created: {normal_output_file}")
else:
    print(f"Failed to create filtered files for {cancer_file}")
    sys.exit(1)



##############################################################################################

# Part 5: Copy files to 'dist' directory and change working directory
dist_dir = 'dist'
if not os.path.exists(dist_dir):
    os.makedirs(dist_dir)

for file in os.listdir('.'):
    if file.startswith('filtered_') and file.endswith('_cancer.parquet'):
        shutil.copy(file, os.path.join(dist_dir, file))
    elif file.endswith('_comb.vcf_procs.parquet'):
        shutil.copy(file, os.path.join(dist_dir, file))

for file in os.listdir('.'):
    if file.startswith('filtered_') and file.endswith('_cancer.parquet'):
        shutil.copy(file, os.path.join(dist_dir, file))
        if os.path.exists(os.path.join(dist_dir, file)):
            print(f"File copied to 'dist': {file}")
        else:
            print(f"Failed to copy file to 'dist': {file}")
            sys.exit(1)

os.chdir(dist_dir)

##################################################################################################

# Part 6: Process files in 'dist' directory

# Check if 'dist' directory contains the required files
if not os.listdir('.'):
    print(colored("No files found in 'dist' directory. Exiting Part 6.", "red"))
    exit()

for file in os.listdir('.'):
    if file.endswith('_comb.vcf_procs.parquet'):
        # Construct the filtered file name based on the original file name
        parts = file.split('_')
        filtered_file = f'filtered_{parts[0]}_{parts[1]}_cancer.parquet'

        if os.path.exists(filtered_file):
            # Load the data
            original_data = pd.read_parquet(file)
            filtered_data = pd.read_parquet(filtered_file)

            # Calculate variant frequencies in the original dataset
            original_data['Variant_Frequency'] = original_data['Individual_1'].apply(lambda x: x.count('1') / len(x))

            # Calculate variant frequencies in the filtered dataset
            filtered_data['Variant_Frequency'] = filtered_data['Individual_1'].apply(lambda x: x.count('1') / len(x))

            # Obtain summary statistics
            original_summary = original_data['Variant_Frequency'].describe()
            filtered_summary = filtered_data['Variant_Frequency'].describe()

            # Print the summaries to the console
            print("Original Data Summary:")
            print(original_summary)
            print("\nFiltered Data Summary:")
            print(filtered_summary)

            # Plotting histograms for visual comparison
            plt.figure(figsize=(12, 6))

            plt.subplot(1, 2, 1)
            plt.hist(original_data['Variant_Frequency'], bins=30, color='blue', alpha=0.7)
            plt.title('Original Variant Frequencies')
            plt.xlabel('Frequency')
            plt.ylabel('Count')

            plt.subplot(1, 2, 2)
            plt.hist(filtered_data['Variant_Frequency'], bins=30, color='green', alpha=0.7)
            plt.title('Filtered Variant Frequencies')
            plt.xlabel('Frequency')
            plt.ylabel('Count')

            # Save the figure
            graph_filename = file.replace('.parquet', '') + '_distribution_comparison.jpeg'
            plt.tight_layout()
            plt.savefig(graph_filename)
            plt.close()

            # Save the distribution summaries to a text file
            summary_filename = file.replace('.parquet', '') + '_distribution_summaries.txt'
            with open(summary_filename, 'w') as file:
                file.write("Original Data Summary:\n")
                file.write(str(original_summary))
                file.write("\n\nFiltered Data Summary:\n")
                file.write(str(filtered_summary))

            print(f"Graph saved as {graph_filename}")
            print(f"Summary saved as {summary_filename}")

print(colored("Analysis and plotting for files in 'dist' directory completed.", "yellow"))

##################################################################################################

# Part 7: Create 'imputer' subdirectory in the parent directory and copy filtered files

# Ensure we are in the parent directory where filtered files are located
parent_directory = os.path.dirname(os.getcwd())
os.chdir(parent_directory)
print(f"Changed to parent directory: {os.getcwd()}")

# Create 'imputer' subdirectory if it doesn't exist
imputer_dir = 'imputer'
if not os.path.exists(imputer_dir):
    os.makedirs(imputer_dir)
    print(f"'imputer' directory created at: {os.path.join(parent_directory, imputer_dir)}")

# Copy all filtered cancer and normal files to 'imputer' directory
for file in glob.glob('filtered_*.parquet'):
    source_path = file
    destination_path = os.path.join(imputer_dir, file)
    shutil.copy(source_path, destination_path)
    if os.path.exists(destination_path):
        print(f"File copied to 'imputer': {file}")
    else:
        print(f"Failed to copy file to 'imputer': {file}")
        sys.exit(1)

print(colored("Filtered and normal files copied to 'imputer' directory.", "yellow"))



#################################################################################################

# Part 8: Apply imputation and split the data into training and testing sets

# Change current working directory to 'imputer' directory
current_directory = os.getcwd()
imputer_directory = os.path.join(current_directory, 'imputer')
os.chdir(imputer_directory)

print(colored("Starting imputation and data splitting.", "magenta"))

# Function to extract chromosome and cancer type from filename
def extract_info_from_filename(filename):
    # Assuming the format is 'filtered_ch{chromosome}_{cancer_type}.parquet'
    parts = filename.split('_')
    chromosome = parts[1]  # Assuming the chromosome is the second element
    cancer_type = parts[2].split('.')[0]  # Assuming the cancer type is the third element and removing file extension
    return chromosome, cancer_type

# Initialize an imputer with mean strategy for numeric data
imputer = SimpleImputer(missing_values=np.nan, strategy='mean')

# Scanning the current directory for all files that begin with 'filtered_' and end with '.parquet'
for file_path in glob.glob('./filtered_*.parquet'):
    print(colored(f"Processing file: {file_path}", 'magenta'))
    data = pd.read_parquet(file_path)

    # Determine numeric columns
    numeric_cols = data.select_dtypes(include=[np.number]).columns

    # Apply imputation to numeric columns only
    data[numeric_cols] = imputer.fit_transform(data[numeric_cols])

    # Splitting data into training and testing sets (70% train, 30% test)
    X_train, X_test = train_test_split(data, test_size=0.3, random_state=42)

    # Extracting chromosome and cancer type from the filename
    chromosome, cancer_type = extract_info_from_filename(os.path.basename(file_path))
    
    # Construct new file names
    train_file_name = f'imputed_{chromosome}_{cancer_type}_train.parquet'
    test_file_name = f'imputed_{chromosome}_{cancer_type}_test.parquet'

    # Save the training and testing data back to new .parquet files
    X_train.to_parquet(train_file_name)
    X_test.to_parquet(test_file_name)

    # Display partial results in white for the training set
    print(colored(f"Partial results for {train_file_name}:\n{X_train.head()}", 'white'))
    print(colored(f"Saved imputed training data to {train_file_name}", 'magenta'))
    print(colored(f"Saved imputed testing data to {test_file_name}\n", 'magenta'))

        
################################################################################################

# Part 9: Analyze the cancer and normal data using Logistic Regression

print(colored("Starting cancer and normal data analysis.", "magenta"))

def main():
    # Function to load data and extract features, labels, and identifiers
    def load_data(file_path, label):
        data = pd.read_parquet(file_path)
        features = data[['Variant_Frequency']].copy()
        identifiers = data[['Chromosome', 'Position']]
        features['label'] = label
        return features, identifiers

    # Function to save results to a text file
    def save_results_to_file(file_name, content):
        with open(file_name, "w") as text_file:
            text_file.write(content)

    # Function to generate file name based on cancer file naming convention
    def generate_file_name(chromosome, cancer_type, file_type):
        if file_type == "model_results":
            return f"model_results_ch{chromosome}_{cancer_type}.txt"
        else:
            return f"cancer_associated_variants_ch{chromosome}_{cancer_type}.csv"

    # Function to extract chromosome number and cancer type from file name
    def extract_info(file_name):
        match = re.search(r'imputed_ch(\d+|X|Y)_(\w+)_', file_name, re.IGNORECASE)
        return (match.group(1).upper(), match.group(2)) if match else (None, None)

    # Automatically find the cancer and normal file paths
    cancer_files = glob.glob('imputed_ch*_*.parquet')
    normal_train_files = glob.glob('imputed_*_normal_train.parquet')
    normal_test_files = glob.glob('imputed_*_normal_test.parquet')

    if not cancer_files:
        print(colored("Cancer files not found. Please check the file naming and location.", "red"))
        return

    # Process each cancer file
    for cancer_file in cancer_files:
        chromosome, cancer_type = extract_info(cancer_file)
        if not chromosome or not cancer_type:
            print(colored(f"Invalid file name format: {cancer_file}", "yellow"))
            continue

        is_train = 'train' in cancer_file
        file_type = 'train' if is_train else 'test'
        cancer_file_path = cancer_file
        normal_file = next((f for f in (normal_train_files if is_train else normal_test_files) if extract_info(f)[0] == chromosome), None)

        if not normal_file:
            print(colored(f"Normal file for chromosome {chromosome} and type {file_type} not found.", "yellow"))
            continue

        # Load the data for cancer and normal cohorts
        print(colored(f"Loading data for chromosome {chromosome}, cancer type {cancer_type}, type {file_type}...", "green"))
        cancer_data, cancer_ids = load_data(cancer_file_path, 1)
        normal_data, normal_ids = load_data(normal_file, 0)

        # Combine cancer and normal data for training and testing
        X_train = pd.concat([cancer_data, normal_data])
        y_train = X_train.pop('label')  # Extracting the labels
        test_identifiers = pd.concat([cancer_ids, normal_ids])  # Combined identifiers for test set

        # Applying SMOTE to the training data
        print("Applying SMOTE to training data...")
        smote = SMOTE()
        X_train_smote, y_train_smote = smote.fit_resample(X_train, y_train)

        # Initialize and train logistic regression model
        model = LogisticRegression()
        model.fit(X_train_smote, y_train_smote)

        # Predictions
        y_pred = model.predict(X_train)

        # Evaluating the model
        accuracy = accuracy_score(y_train, y_pred)
        class_report = classification_report(y_train, y_pred)

        # Print results to the screen
        print("Accuracy:", accuracy)
        print("Classification Report:")
        print(class_report)

        # Save results to a text file
        results_text = f"Accuracy: {accuracy}\nClassification Report:\n{class_report}"
        results_file_name = generate_file_name(chromosome, cancer_type, "model_results")
        save_results_to_file(results_file_name, results_text)

        print(f"Results have been saved to {results_file_name}")

        # Combine predictions with test dataset identifiers to list variants
        predicted_data = test_identifiers.copy()
        predicted_data['Predicted_Label'] = y_pred  # Add predictions as a new column

        # Extracting cancer-associated variants
        cancer_associated_variants = predicted_data[predicted_data['Predicted_Label'] == 1]

        # Save the cancer-associated variants to a file including 'Chromosome' and 'Position'
        variants_file_name = generate_file_name(chromosome, cancer_type, "variants")
        cancer_associated_variants.to_csv(variants_file_name, index=False)

        print(f"Cancer-associated variants with chromosome and position have been saved to {variants_file_name}")

      

if __name__ == "__main__":
    main()

    
print(colored("Cancer and normal data analysis completed.", "yellow"))




################################################################################################

# Part 10: Create 'snps' directory in the parent directory and copy specific files to it

# Get the parent directory of the current working directory
parent_directory = os.path.dirname(os.getcwd())
print(f"Parent directory: {parent_directory}")

# Create the 'snps' directory in the parent directory
snps_dir = os.path.join(parent_directory, 'snps')
if not os.path.exists(snps_dir):
    os.makedirs(snps_dir)
    print(f"'snps' directory created at: {snps_dir}")

# Copy all 'cancer_associated_variants_' files to the 'snps' directory
for file in glob.glob(os.path.join(os.getcwd(), 'cancer_associated_variants_*.csv')):
    source_path = file
    destination_path = os.path.join(snps_dir, os.path.basename(file))
    shutil.copy(source_path, destination_path)
    if os.path.exists(destination_path):
        print(f"File copied to 'snps': {os.path.basename(file)}")
    else:
        print(f"Failed to copy file to 'snps': {os.path.basename(file)}")
        sys.exit(1)

print(colored("All specified files copied to 'snps' directory.", "yellow"))


################################################################################################

# Part 11: Obtain SNP accession numbers and save to new file


# Function to save progress
def save_progress(df, directory, filename, last_index):
    full_path = os.path.join(directory, filename)
    df.to_csv(full_path, index=False)
    with open(f"{full_path}_last_index.txt", 'w') as f:
        f.write(str(last_index))

# Function to get the last processed index
def get_last_processed_index(filename):
    try:
        with open(f"{filename}_last_index.txt", 'r') as f:
            return int(f.read().strip())
    except FileNotFoundError:
        return -1

# Function to send email via SendGrid
def send_email_via_sendgrid(from_email, to_email, subject, content):
    api_key = os.getenv('SENDGRID_API_KEY')
    if not api_key:
        print("SendGrid API key is not set. Check your environment variables.")
        return

    sg = sendgrid.SendGridAPIClient(api_key=api_key)
    from_email = Email(from_email)
    to_email = To(to_email)
    content = Content("text/plain", content)

    mail = Mail(from_email, to_email, subject, content)
    try:
        response = sg.client.mail.send.post(request_body=mail.get())
        print("Email sent! Status Code:", response.status_code)
    except Exception as e:
        print(f"Failed to send email: {e}")

def main():
    # Define the base directory relative to this script's location
    base_directory = os.path.dirname(os.path.realpath(__file__))
    snps_dir = os.path.join(base_directory, 'snps')

    if not os.path.exists(snps_dir):
        print(f"Directory '{snps_dir}' not found.")
        sys.exit(1)

    for file in glob.glob(os.path.join(snps_dir, 'cancer_associated_variants_ch*.csv')):
        match = re.search(r'cancer_associated_variants_ch(\d+)_([a-zA-Z0-9]+)', os.path.basename(file))
        if not match:
            print(f"File name format is incorrect: {file}")
            continue

        chromosome_number, cancer_type = match.groups()
        file_path = os.path.join(snps_dir, os.path.basename(file))
        job_title = f"Chromosome {chromosome_number} {cancer_type}"
        user_email = "user@example.com"  # Replace with your email

        print(f"Analyzing file: {os.path.basename(file_path)}")

        df = pd.read_csv(file_path)
        if df.empty:
            print(f"The DataFrame for {file_path} is empty. Skipping.")
            continue

        df.dropna(subset=['Chromosome', 'Position'], inplace=True)
        df['Chromosome'] = df['Chromosome'].astype(int)
        df['Position'] = df['Position'].astype(int)

        if 'SNP' not in df.columns:
            df['SNP'] = 'N/A'

        new_file_name = f"snp_access_output_chr{chromosome_number}_{cancer_type}.csv"
        last_processed_index = get_last_processed_index(os.path.join(snps_dir, new_file_name))

        total_rows = len(df)
        rows_left = total_rows - last_processed_index - 1
        print(f"Total rows to process: {rows_left}")

        print('\033[95mFetching SNP accession numbers from Ensembl...\033[0m')
        
        try:
            for index, row in df.loc[last_processed_index+1:].iterrows():
                rows_left = total_rows - index - 1
                print(f"Processing row {index + 1}/{total_rows} (Rows left: {rows_left})")

                endpoint = f"http://rest.ensembl.org/overlap/region/human/{row['Chromosome']}:{row['Position']}-{row['Position']}?feature=variation"
                response = requests.get(endpoint, headers={"Content-Type": "application/json"})

                if response.status_code == 200 and response.json():
                    snp_id = response.json()[0]['id']
                    df.at[index, 'SNP'] = snp_id
                    print(f"Position: {row['Position']}, SNP: {snp_id}")
                else:
                    print(f"Position: {row['Position']}, SNP: N/A")

                if index % 100 == 0:
                    save_progress(df, snps_dir, new_file_name, index)

            save_progress(df, snps_dir, new_file_name, df.index[-1])
            content = "The analysis has completed successfully."
            send_email_via_sendgrid('sudoroot1775@outlook.com', user_email, f"{job_title}_Analysis Completed", content)
            print(f'\033[95mSNP accession numbers saved as {new_file_name}.\033[0m')

        except Exception as e:
            print(f"An error occurred: {e}")
            content = f"An error occurred during the analysis: {str(e)}"
            try:
                send_email_via_sendgrid('sudoroot1775@outlook.com', user_email, f"{job_title}_Analysis Interrupted", content)
            except Exception as email_error:
                print(f"Failed to send interruption email: {email_error}")
            # Do not exit, continue with the next file
            print(f"Continuing with the next file.")

if __name__ == "__main__":
    main()


