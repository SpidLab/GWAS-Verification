# GWAS-Verification
The source code is an artifact of the paper "Privacy-Preserving and Efficient Verification of the Outcome in Genome-Wide Association Studies" which can be accessed in https://arxiv.org/abs/2101.08879. 
At the time of deployment, we used Python version '3.7.9'(for the main implementation) and Matlab '9.11.0.1769968 (R2021b)' (to find the cutoff points) along with some other packages which are included in *requirements.txt*.

In the paper we used different phenotypes of openSNP data, which are not provided due to privacy policies. Nevertheless, they can be crawled from https://opensnp.org/phenotypes. After gathering the data, one may use *preprocess_dataset.py* to preprocess the data and convert it to the file format used in the paper (check "Data/D_original_dataset.csv" as a reference format). For the next steps, we use toy data to show the computations performed bye the two main parties: **researcher** and **verifier**. The parameters and other details to run each of the scripts are explained below.

## preprocess_dataset.py
The main two inputs of this script include the raw data from each sample (user) and the user IDs. Taking in consideration the huge number of SNPs per user, generally only the data of one chromosome is extracted. Thus, it is necessary to provide the chromosome number (1-22). And, since there might be different SNPs in different users' data for a particular chromosome (due to sequencing errors), it is needed to determine which SNPs of that chromosome will be used. Therefore, we also provide one the user's file in order to extract the SNPs. Lastly, there is a parameter to specify the output directory for the constructed dataset. 

__Required parameters__

Flag | Description 
--- | ---
-i data_dir |  Path to the directory which contains the users (samples)
-u user_IDs_file | A single line comma separated csv file that contains user IDs
-o output_dir | Path where the output file will be saved
-c chromosome_number | Determines the data of which chromosome will be extracted
-s user_SNPs_file | Full path to the user from which the SNP IDs will be selected


Example run
```
python preprocess_dataset.py -i "Data/Users/" -u "Data/user_IDs.csv" -o "Data/" -c 1 -s "Data/Users/1"
```

## researcher_computations.py
This script illustrates the functionalities of a researcher. After taking as an input the dataset *D* and the user IDs for case and control groups, the researcher initially performs GWAS on dataset *D*, sorts them ascendingly according to the p-value and later adds noise to top *k* SNPs. In the code, we assume that the researcher provides the SNP IDs with indexes from *0* to *l*, but shifts down on the SNP list with a parameter *shift* and provides the partial noisy dataset of the SNPs with indexes from *shift* to *shift + l* (corresponding to underselling scenario in the paper). Note that when *shift = 0*, the researcher provides the correct statistics. One may change these indexes to try overselling and mixed scenarios mentioned in the paper). Finally, it outputs the GWAS results for *l* SNPs (R<sup>t</sup><sub>l</sub>) and the partial noisy dataset (D<sup>ε</sup><sub>k</sub>) in a pre-specified folder.


__Required parameters__

Flag | Description 
--- | ---
-i dataset_file |  Path to the original dataset file *D*
-e epsilon | Privacy parameter - ε
-k k | Number of SNPs that are provided in the partial noisy dataset
-l l | Number of SNPs for which GWAS statistics are provided
-s shift | Determines the starting point for which the SNPs statistics are provided
-c D_case_control_IDs_file | Path to the csv file that contains case and control IDs for dataset *D*
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
