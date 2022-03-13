# GWAS-Verification
Description of the project

## generate_dataset.py

__Required parameters__

Flag | Description 
--- | ---
-i DATA_PATH |  Path to the directory which contains the users (samples).
-c CHROMOSOME_NUMBER | Pick for which chromosome will be extracted the data.
-s SEQUENCED_USER_ID | Full path to the user from which the SNP IDs will be selected.
-u USER_IDS_FILE | Path to the csv files that contains case and control IDs.
-o OUTPUT_FILE | Path to save the output file.

See the --help
```
usage: calculate_kinship.py [-h] -i DATA_PATH -c CHROMOSOME_NUMBER -s SEQUENCED_USER_ID -u USER_IDS_FILE -o OUTPUT_FILE

Generate dataset

optional arguments:
  -h, --help            show this help message and exit
  -i DATA_PATH, --data_path DATA_PATH
                        Path to the directory which contains the users
                        (samples).
  -c CHROMOSOME_NUMBER, --chromosome_number CHROMOSOME_NUMBER
                        Pick for which chromosome will be extracted the data.
  -s SEQUENCED_USER_ID, --sequenced_user_ID SEQUENCED_USER_ID
                        Path to the user from which the SNP IDs will be
                        selected.
  -u USER_IDS_FILE, --user_IDs_file USER_IDS_FILE
                        Path to the csv files that contains case and control
                        IDs.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to save the output file.
```

Example run
```
>python generate_dataset.py -i "Data/Users/" -c "22" -s 11 -u "Data/userIDs.csv" -o "Data/"
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
