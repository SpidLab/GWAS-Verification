import csv
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Generate dataset")
parser.add_argument('-i', '--data_dir', type=str, help='Path to the directory which contains the users (samples).', required=True)
parser.add_argument('-u', '--user_IDs_file', type=str, help='A single line comma separated csv file that contains user IDs', required=True)
parser.add_argument('-o', '--output_dir', type=str, help='Path where the output file will be saved', required=True)
parser.add_argument('-c', '--chromosome_number', type=str, help='Determines the data of which chromosome will be extracted', required=True)
parser.add_argument('-s', '--user_SNPs_file', type=str, help='Full path to the user from which the SNP IDs will be selected', required=True)
args = parser.parse_args()

def extract_SNPs():
    SNP_list = [] # holds the list of SNPs
    with open(args.user_SNPs_file, "r", encoding="utf8") as tsv:
        for line in csv.reader(tsv, delimiter="\t"):
            if not (line[0].startswith("#")) and len(line) > 2:
                if args.chromosome_number == "NA":
                    SNP_list.append(line[0])
                else:
                    if str(line[1]) == str(args.chromosome_number):
                        SNP_list.append(line[0])
    return SNP_list

def get_user_IDs():
    with open(args.user_IDs_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for index, row in enumerate(spamreader):
            user_IDs = row
    return user_IDs

def construct_dataset_from_files(SNP_list, user_IDs):
    file_index = 0  # holds the indexing of which file is being processed
    inverted_dataset = []  # holds a 2D matrix where a row represents a user and a column represents a SNP (the inverted representation)

    for user_ID in user_IDs:  # iterate over each user ID
        print(str(file_index+1) + " out of " + str(len(user_IDs)) + " users are processed")
        full_path = args.data_dir + user_ID

        dictionary = {} # holds the allele for each SNP (e.g. dictionary = {'rs3905942': 'GG', ....}
        with open(full_path, "r", encoding="utf8") as tsv:
            for line in csv.reader(tsv, delimiter="\t"):
                if len(line) > 2 and not (line[0].startswith("#")) and not (line[0].startswith("rsid")): # avoids some common errors in files
                    if args.chromosome_number == "NA":
                        dictionary[line[0]] = line[3]
                    else:
                        if str(line[1]) == str(args.chromosome_number):
                            dictionary[line[0]] = line[3]

        # Generate an empty 2d matrix
        row = []
        for x in range(len(SNP_list)):
            row.append("-")
        inverted_dataset.append(row)

        # Fill the values in the 2d matrix
        for SNP_index, value in enumerate(SNP_list, start=0):
            if value in dictionary:
                inverted_dataset[file_index][SNP_index] = dictionary[value]
        file_index += 1
    return inverted_dataset

def convert_df_to_numbers(genomic_df):
    genomic_df_values = genomic_df.values # gets the values inside the dataframe
    major_SNPs = [] #array that holds the major allele for each SNP
    minor_SNPs = [] #array that holds the minor allele for each SNP

    for i, row in enumerate(genomic_df_values, start=0):
        freq = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        } # holds the frequency of each nucleotide per SNP
        for j, entry in enumerate(row, start=0):
            for char in entry:
                if (char == "A"):
                    freq["A"] += 1
                elif (char == "C"):
                    freq["C"] += 1
                elif (char == "G"):
                    freq["G"] += 1
                elif (char == "T"):
                    freq["T"] += 1

        sorted_freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)
        major = sorted_freq[0][0]  # holds the most frequent nucleotide
        minor = sorted_freq[1][0]  # 2nd most frequent nucleotide

        for j, entry in enumerate(row, start=0):
            minor_count = 0  # can be 0, 1 or 2 dependant on nr of minor allele
            for char in entry:
                if (char == minor):
                    minor_count += 1
            # modify the dataframe from strings (e.g AG) to numbers (0, 1, 2)
            genomic_df_values[i][j] = minor_count
        major_SNPs.append(major)
        minor_SNPs.append(minor)
    # Convert the values to the numbered pandas dataframe
    numbered_df = pd.DataFrame(data=genomic_df_values, index=genomic_df.index, columns=genomic_df.columns)

    numbered_df['alt'] = minor_SNPs #a column holding alternate (minor) alleles
    numbered_df['ref'] = major_SNPs #a column holding reference (major) alleles
    return  numbered_df

if __name__ == "__main__":
    # Selects the SNPs for which the whole dataset will be built from one of the sequenced samples
    SNP_list = extract_SNPs()
    print("The total number of SNPs is: " + str(len(SNP_list)))

    # Get user IDs
    user_IDs = get_user_IDs()
    print("The total number of samples is: " + str(len(user_IDs)))

    # Generate the dataset (which is initially inverted) from the user files (samples)
    inverted_dataset = construct_dataset_from_files(SNP_list, user_IDs)

    # Generate pandas dataframe from the 2d matrix and add SNP_list to it, then transpose the matrix
    print("Converting the genomic dataframe to a numbered dataframe")
    genomic_df = pd.DataFrame(inverted_dataset, columns=SNP_list)
    genomic_df.index = user_IDs
    genomic_df = genomic_df.T #transpose

    numbered_df = convert_df_to_numbers(genomic_df)
    numbered_df.to_csv(args.output_dir + "original_dataset.csv")
    print(str(args.output_dir) + "original_dataset.csv generated successfully")