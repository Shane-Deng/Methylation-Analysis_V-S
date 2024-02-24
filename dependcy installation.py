import subprocess

env_name = str(input("What do you want to call your environment").strip())
# Create and activate a new Conda environment (Modify 'myenv' and Python version as needed)
python_version = str(input("Which python version do you want? 3.8 is recommended").strip())
subprocess.run(f"conda create --name {env_name} python={python_version} -y", shell=True)
subprocess.run(f"conda activate {env_name}", shell=True)

# Install Python packages with conda install
python_packages = [
    "pyfiglet",
    "termcolor",
    "sendgrid"
]
for package in python_packages:
    subprocess.run(f"conda install -c conda-forge {package} -y", shell=True)

# Install external bioinformatics tools with conda install
bioinformatics_tools = [
    "bismark",
    "bowtie2",
    "samtools"
]
for tool in bioinformatics_tools:
    subprocess.run(f"conda install -c bioconda {tool} -y", shell=True)

print(f"Installation completed. Please do the following step:\n 1. run "conda init bash" \n2. run"conda activate {env_name}")


