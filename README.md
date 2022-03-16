# GWAS-Verification
The source code is an artifact of the paper "Privacy-Preserving and Efficient Verification of the Outcome in Genome-Wide Association Studies" which can be accessed in https://arxiv.org/abs/2101.08879. 
At the time of deployment, we used Python version '3.7.9'(for the main implementation) and Matlab '9.11.0.1769968 (R2021b)' (to find the cutoff points) along with corresponding packages which are included in *requirements.txt*.

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
This script illustrates the functionalities of a researcher. After taking as an input the dataset *D* and the user IDs for case and control groups, the researcher initially performs GWAS on dataset *D*, sorts them ascendingly according to the p-value and later adds noise to top *k* SNPs using a provided ε value. In the code, we assume that the researcher provides the SNP IDs with indexes from *0* to *l*, but shifts down on the SNP list with a parameter *shift* and provides the partial noisy dataset of the SNPs with indexes from *shift* to *shift + l* (corresponding to underselling scenario in the paper). Note that when *shift = 0*, the researcher provides the correct statistics. One may change these indexes to try overselling and mixed scenarios mentioned in the paper). Finally, it outputs the GWAS results for *l* SNPs (R<sup>t</sup><sub>l</sub>) and the partial noisy dataset (D<sup>ε</sup><sub>k</sub>) in a pre-specified folder.

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

This script illustrates the functionalities of a verifier.
Verifier initially takes as input the dataset *E*, the corresponding user IDs for case and control groups, GWAS statistics of dataset *D* and partial noisy dataset (R<sup>t</sup><sub>l</sub> and D<sup>ε</sup><sub>k</sub> which were previously generated by researcher), and the user IDs for case and control groups of dataset *D*. Furthermore, it also requires as inputs the privacy parameter ε, and pre-computed threshold values which are needed to classify SNPs as correct or incorrect based on odds ratio, MAF and p-value (τ<sub>o</sub>, τ<sub>p</sub>,  τ<sub>a</sub>). After doing the necessary computations (such as finding the deviation and the error between results), the verifier outputs a dataframe (a table alike data structure in Python), where SNPs are represented by rows and for each SNP it outputs the error and correctness classification for each statistic (odds ratio, MAF, p-value).

__Required parameters__

Flag | Description 
--- | --- 
-i dataset_file |  Path to the original dataset file (E)
-a E_case_control_IDs_file |  Path to the csv file that contains case and control IDs for dataset E
-n noisy_dataset_file |  Path to the partial noisy dataset D<sup>ε</sup><sub>k</sub>
-g D_dataset_GWAS_file |  Path to the GWAS of dataset *D* file
-b D_case_control_IDs_file |  Path to the csv file that contains case and control IDs for dataset D
-e epsilon | Privacy parameter - ε
-x odds_cut_off | Threshold (cut-off point) for odds ratio - τ<sub>o</sub>
-y maf_cut_off | Threshold (cut-off point) for MAF - τ<sub>p</sub>
-z pval_cut_off | Threshold (cut-off point) for p-value - τ<sub>a</sub>
-o output_dir | Directory to save the output files

Example run
```
python verifier_computations.py -i "Data/E_original_dataset.csv" -a "Data/E_case_control_IDs.csv" -n "Data/Researcher/noisy_dataset_D.csv" -g "Data/Researcher/D_GWAS.csv" -b "Data/D_case_control_IDs.csv" -e 3 -x 0.8 -y 0.00001 -z 0.5 -o "Data/Verifier/"
```

## cutoff_computation.py

In order to compute the cutoff point, the verifier needs as input datasets *E* and *F*, along with the user IDs for their case and control groups, ε value, *l*, *shift* (which are all explained above) and a number of splits which determines how many different values will be on the x-axis for the probability distribution plotting. If sufficient data is generated, increasing nr_splits generally results in finding more accurate cutoff points. To generate the line for correct statistics, value of *shift* needs to be set to 0. 
```
python cutoff_computation.py -a "Data/E_case_control_IDs.csv" -b "Data/F_case_control_IDs.csv" -c "Data/E_original_dataset.csv" -d "Data/F_original_dataset.csv" -e 3 -s 0 -l 100 -n 10
```
For incorrect statistics, picking different *shift* values results in different deviation of the provided statistic with the correct ones. Below is an example run for *shift* = 100: 
```
python cutoff_computation.py -a "Data/E_case_control_IDs.csv" -b "Data/F_case_control_IDs.csv" -c "Data/E_original_dataset.csv" -d "Data/F_original_dataset.csv" -e 3 -s 100 -l 100 -n 10
```
For a more accurate cut-off point, a verifier may also combine the results of picking different *shift* values (and normalizing the combined results accordingly) before plotting them on MATLAB. The two above commands generate "index.txt", "correct.txt" and "incorrect.txt", which are sent to MATLAB folder to be plotted. In *cutoff_computation.py* script, we only show how to find the cutoff point for p-value, but it can be easily modified to compute the cutoff value for other statistics as well.

__Required parameters__

Flag | Description 
--- | --- 
-a E_case_control_IDs_file |  Path to the csv file that contains case and control IDs for dataset *E*
-b F_case_control_IDs_file |  Path to the csv file that contains case and control IDs for dataset *F*
-c E_dataset_file |  Path to the csv file that contains dataset *E*
-d F_dataset_file |  Path to the csv file that contains dataset *F*
-e epsilon | Privacy parameter - ε
-s shift | Determines the starting point for which the SNPs statistics are provided
-l l | Number of SNPs for which GWAS statistics are provided
-n nr_splits | Plot parameter to define number of split on the x-axis

## intersection.m

This scripts plots the two lines (correct and incorrect statistics) from the provided files and finds the intersection point between them which corresponds to the cutoff point for p-value. In the current version, a toy example is provided as the input to the script to show what a possible plot would look like (the toy example filenames needs to be replaced with the commented out lines).
