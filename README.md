# GWAS-Verification
Description of the project

## preprocess_dataset.py

__Required parameters__

Flag | Description 
--- | ---
-i data_dir |  Path to the directory which contains the users (samples).
-u user_IDs_file | Path to the csv files that contains user IDs.
-o output_dir | Path to save the output file.
-c chromosome_number | Pick for which chromosome will be extracted the data.
-s user_SNPs_file | Full path to the user from which the SNP IDs will be selected.


Example run
```
>python generate_dataset.py -i "Data/Users/" -u "Data/user_IDs.csv" -o "Data/" -c "22" -s "Data/Users/11"
```

## researcher_operations.py

__Required parameters__

Flag | Description 
--- | ---
-p DATASET_PATH |  Path to the directory which contains the users (samples).
-e EPSILON | Epsilon value
-k K | Number of SNPs that are provided in the partial noisy dataset
-l L | Number of SNPs that are provided as part of GWAS statistics
-f FLAG | Whether to report correct or incorrect statistics
-u USER_IDS_FILE | Path to the csv files that contains case and control IDs.
-o OUTPUT_DIR | Directory to save the output files
