# Full-ConDoR Pipeline

This repository contains a workflow for running the Full-ConDoR pipline using Snakemake. 

## Prerequisites

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Bolladeen/full-ConDoR.git
   cd full-ConDoR
   ```

2. **To run the pipeline, create a Conda environment with the required dependencies using the command below:**

    ```bash
    conda create -n condor_env python=3.8 snakemake pandas numpy matplotlib seaborn pyyaml scipy biopython ete3 scikit-learn click networkx gurobipy
    conda activate condor_env
    ```

3. **Generate Config File**

    ```yaml
    # Example configuration file
    
    # Define the working directory for output files
    workdir: /path/to/working/directory
    
    # List of patients to process
    patients:
      - patient1
      - patient2
    
    # Path to raw data directory
    raw_data_directory: /path/to/raw_data
    
    # Path to the Falcon solutions directory
    falcon_solutions: /path/to/falcon_solutions
    
    # Path to the directory with annotated mutations
    annotated_mutations: /path/to/annotated_mutations
    
    # Path to the directory with subclonal mutations (optional)
    subclonal_mutations: /path/to/subclonal_mutations
    
    # Path to the Condor pipeline scripts directory
    condor_pipeline_scripts_dir: /path/to/condor_pipeline_scripts
    
    # Path to the amplicon coordinates file
    amplicon_coordinates_file: /path/to/amplicon_coordinates.txt
    
    # Path to the downstream Condor scripts directory
    condor_downstream_scripts_dir: /path/to/condor_downstream_scripts
    
    # Path to the fast Condor script
    fast_condor_script: /path/to/fast_condor_script.py
    
    # Path to the sample name map file
    sample_name_map: /path/to/sample_name_map.txt
    
    ```
An example config file can be found at ```config_MSK_hz_local.yaml```. 

#### Parameters

- **`workdir`**: Defines the directory where all output files will be generated.

- **`patients`**: A list of patient IDs to process.

- **`raw_data_directory`**: Location of raw patient data in `.h5` format.

- **`falcon_solutions`**: Directory containing Falcon solutions such as clone assignments and profiles.

- **`annotated_mutations`**: Directory of mutation files for each patient.

- **`subclonal_mutations`**: (Optional) Directory for subclonal mutations in YAML format.

- **`condor_pipeline_scripts_dir`**: Directory containing pipeline scripts for Condor.

- **`amplicon_coordinates_file`**: File mapping gene amplicons to their coordinates.

- **`condor_downstream_scripts_dir`**: Directory for downstream Condor scripts (tree generation, etc.).

- **`fast_condor_script`**: Script to run the fast Condor pipeline.

- **`sample_name_map`**: Maps sample names to their corresponding patient IDs.

  

4. **Run Snakemake**
Navigate to the directory containing the condor_pipeline.smk file and run Snakemake using the command:
    ```bash
    snakemake -s condor_pipeline.smk
    ```
