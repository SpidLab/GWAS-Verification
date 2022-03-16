# GWAS-Verification
The source code is an artifact of the paper "Privacy-Preserving and Efficient Verification of the Outcome in Genome-Wide Association Studies" which can be accessed in https://arxiv.org/abs/2101.08879. 

In the paper we used different phenotypes of openSNP data, which are not provided due to privacy policies. Nevertheless, they can be crawled from https://opensnp.org/phenotypes. After gathering the data, one may use *preprocess_dataset.py* to preprocess the data and convert it to the file format used in the paper (check "Data/D_original_dataset.csv" as a reference format). For the next steps, we use toy data to show the computations performed bye the two main parties: **researcher** and **verifier**.


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
python preprocess_dataset.py -i "Data/Users/" -u "Data/user_IDs.csv" -o "Data/" -c "22" -s "Data/Users/11"
```

## researcher_computations.py

__Required parameters__

Flag | Description 
--- | ---
-i dataset_file |  Path to the original dataset file (D) | D<sup>ε</sup>
-e epsilon | Privacy parameter
-k k | Number of SNPs that are provided in the partial noisy dataset
-l l | Number of SNPs for which GWAS statistics are provided
-s shift | The index of SNP that is considered the most associated one. If s = 0, the correct SNPs are provided
-c D_case_control_IDs_file | Path to the csv file that contains case and control IDs for dataset D
-o output_dir | Directory to save the output files

Example run
```
python researcher_computations.py -i "Data/D_original_dataset.csv" -e 3 -k 150 -l 100 -s 0 -c "Data/D_case_control_IDs.csv"  -o "Data/Researcher/"
```

## verifier_computations.py

__Required parameters__

Flag | Description 
--- | --- 
-i dataset_file |  Path to the original dataset file (E)
-a E_case_control_IDs_file |  Path to the csv file that contains case and control IDs for dataset E
-n noisy_dataset_file |  Path to the partial noisy dataset (D_k^e)
-g D_dataset_GWAS_file |  Path to the GWAS of dataset D file
-b D_case_control_IDs_file |  Path to the csv file that contains case and control IDs for dataset D
-e epsilon | Privacy parameter
-x odds_cut_off | Pre-defined threshold value to classify SNPs as correct or incorrect based on odds ratio
-y maf_cut_off | Pre-defined threshold value to classify SNPs as correct or incorrect based on MAF
-z pval_cut_off | Pre-defined threshold value to classify SNPs as correct or incorrect based on p-value
-o output_dir | Directory to save the output files



Example run
```
python verifier_computations.py -i "Data/E_original_dataset.csv" -a "Data/E_case_control_IDs.csv" -n "Data/Researcher/noisy_dataset_D.csv" -g "Data/Researcher/D_GWAS.csv" -b "Data/D_case_control_IDs.csv" -e 3 -x 0.8 -y 0.00001 -z 0.5 -o "Data/Verifier/"
```


**Notes:Add which python version I'm using, along with the version package in requirements.txt file**
