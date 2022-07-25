import pandas as pd
import scipy.stats as stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", help="input target genes", required=True)
parser.add_argument("-b", "--background", action="store", help="input background file", required=True)
parser.add_argument("-o", "--output", action="store", help="output file name", required=True)
args = parser.parse_args()

"""
# This script uses Fisher's exact test to determine if a particular chromosome is overrepresented with genes from a gene set of interest
# input data
## input target gene file (no header)
### format: plain text file of gene names
## background file - all genes in genome (no header)
### format: gene_name\tchromosome
## output file
### format: chromosome\tpvalue
"""
# read in target gene dataframe
target_gene_df = pd.read_csv(args.input, sep = '\t', names = ['gene'])

# create a list of target genes 
target_gene_df_list = target_gene_df['gene'].tolist()

# read in background gene dataframe
background_df = pd.read_csv(args.background, sep = '\t', header = None, names = ['gene', 'chromosome'])

# subset target genes from gene chromosome df
target_gene_chromosome_df = background_df[background_df['gene'].isin(target_gene_df_list)]

# initialize empty dataframe list
combined_df_list = []

# iterate over the chromosomes with target genes
for chromosome in sorted(set(target_gene_chromosome_df['chromosome'].tolist())):

    # subset genes that are on chromosome but not target genes
    on_chromosome_not_target = len(background_df[(background_df['chromosome'] == chromosome) & (~background_df['gene'].isin(target_gene_df_list))])

    # subset genes that are on chromosome and are target genes
    on_chromosome_target = len(background_df[(background_df['chromosome'] == chromosome) & (background_df['gene'].isin(target_gene_df_list))])

    # subset genes that are not on chromosome and not target genes
    subset_inv_chromosome_not_target = len(background_df[(background_df['chromosome'] != chromosome) & (~background_df['gene'].isin(target_gene_df_list))])

    # subset genes that are not on chromosome and are target genes
    subset_inv_chromosome_and_target = len(background_df[(background_df['chromosome'] != chromosome) & (background_df['gene'].isin(target_gene_df_list))])

    # fisher's exact test
    oddsratio, pvalue = stats.fisher_exact([[on_chromosome_target, subset_inv_chromosome_and_target],[on_chromosome_not_target, subset_inv_chromosome_not_target]], alternative = 'greater')

    # create data frame from results
    df = pd.DataFrame({'chromosome': [chromosome], 'pvalue': [pvalue]})

    # append dataframe to combined dataframe list
    combined_df_list.append(df)

# concatenate dataframes together
combined_df_fin = pd.concat(combined_df_list, ignore_index=True)

# write out output file
combined_df_fin.to_csv(args.output, index=False, sep = '\t')
