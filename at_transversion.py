from Bio import SeqIO
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", help="input alignment file", required=True)
parser.add_argument("-o", "--output", action="store", help="output tab-delimited AT transversion file", required=True)
args = parser.parse_args()

"""
# This script reads in a multiple sequence alignment of two sequences in fasta format and calculates the AT transversion ratio
## input data
### format: multifasta alignment of two DNA sequences
## output data
### format: seq1_name\tseq2_name\tat_transversion_ratio
"""

# read in multiple sequence alignment
fasta_records = list(SeqIO.parse(args.input, "fasta"))

# store the sequences as variables, convert sequences to uppercase
seq_one = fasta_records[0].seq.upper()
seq_two = fasta_records[1].seq.upper()

# create lists from those sequences
seq_one_list = list(seq_one)
seq_two_list = list(seq_two)

# create a list of AT transversions
at_transversionslist = ["AT", "TA"]

# create an empty sting used to compare nucleotides
compare_seq_string = ""

# create count variables for instances of transversion and missing data ("-")
transversion = 0
count_missing = 0

# for loop to iterate over nucleotides
for seq in zip(seq_one_list,seq_two_list):

    compare_seq_string += seq[0]

    compare_seq_string += seq[1]

    if compare_seq_string in at_transversionslist:
        transversion += 1

    if "-" in compare_seq_string:
        count_missing +=1

    compare_seq_string = ""

# calculate length of alignment minus missing data ("-")
seq_len = len(seq_one) - count_missing

# calculate the AT transversion ratio
at_transversion_ratio = transversion / seq_len

# store sequence names as variables
seq1_name = fasta_records[0].id
seq2_name = fasta_records[1].id

# dictionary which will be converted to df
at_transversion_dict = {'seq1_name':seq1_name,'seq2_name':seq2_name,'at_transversion_ratio': at_transversion_ratio}

# final df with sequence names and AT transversion ratio
at_transversion_df = pd.DataFrame([at_transversion_dict])

# write out file
at_transversion_df.to_csv(args.output, sep='\t', index=False)
