from Bio import AlignIO
import argparse
import re



"""
# This script trims a multifasta alignment to the boundaries of the first sequence in the fasta file
# input
## input multiple sequence alignment
### format: multifasta alignment
# output data
## output multifasta alignment trimmed to the boundaries of the first sequence
### format: multifasta alignment
"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", help="input alignment fasta file", required=True)
parser.add_argument("-o", "--output", action="store", help="output trimmed alignment file", required=True)
args = parser.parse_args()

# read in alignment
aln = AlignIO.read(args.input, 'fasta')

# iterate over sequences and names
seq_list = []
name_list = []
for fasta in aln:
	name, sequence = str(fasta.id), str(fasta.seq)
	seq_list.append(sequence)
	name_list.append(name)

# store first sequence as a variable
first_seq = seq_list[0]

# get index of first letter in first seq
first_seq_letter_index = re.search('[a-z|,]+', first_seq).start()

# get index of last letter in first seq
last_seq_letter_index = list(re.finditer(r'[a-z]', first_seq, re.I))[-1].start() + 1

# trim all sequences using the indexes of the start and end of the first sequence
trim_list = []
for seq in seq_list:
	trim_to_first_seq = seq[first_seq_letter_index:last_seq_letter_index]
	trim_list.append(trim_to_first_seq)

# write out file
i = 0
for a, b in zip(name_list, trim_list):
	name_seq = ">%s\n%s\n" %(a, b)
	if i == 0:
		with open(args.output, "w") as f:
			f.write(name_seq)
			i += 1
	else:
		with open(args.output, "a") as f:
			f.write(name_seq)

	


