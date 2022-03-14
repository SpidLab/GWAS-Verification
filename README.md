# GWAS-Verification
Description of the project

## preprocess_dataset.py

__Required parameters__

Flag | Description 
--- | ---
-i data_dir |  Path to the directory which contains the users (samples)
-u user_IDs_file | Path to the csv files that contains user IDs
-o output_dir | Path to save the output file
-c chromosome_number | Pick for which chromosome will be extracted the data
-s user_SNPs_file | Full path to the user from which the SNP IDs will be selected


Example run
```
>python preprocess_dataset.py -i "Data/Users/" -u "Data/user_IDs.csv" -o "Data/" -c "22" -s "Data/Users/11"
```

## researcher_computations.py

__Required parameters__

Flag | Description 
--- | ---
-i dataset_dir |  Path to the original dataset file
-e epsilon | Privacy parameter
-k k | Number of SNPs that are provided in the partial noisy dataset
-l l | Number of SNPs for which GWAS statistics are provided
-s shift | The index of SNP that is considered the most associated one. If s = 0, the correct SNPs are provided
-c case_control_IDs_file | Path to the csv files that contains case and control IDs
-o output_dir | Directory to save the output files

Example run
```
>python researcher_computations.py -i "Data/original_dataset.csv" -e 3 -k 150 -l 100 -s 0 -c "Data/case_control_IDs.csv"  -o "Data/Researcher/"
```


**
Notes:
Add which python version I'm using, along with the version package in requirements.txt file
Put opensnp link somewhere (https://opensnp.org/phenotypes). **
