import pandas as pd
import pybedtools
from pybedtools import BedTool
pd.options.mode.chained_assignment = None
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", help="input gwas output file", required=True)
parser.add_argument("-p", "--pvalue", action="store", type= float, help="pvalue threshold", required=True)
parser.add_argument("-a", "--annotation_file", action="store", help="gene annotation file", required=True)
parser.add_argument("-b", "--bed_file", action="store", help="gene bed file", required=True)
parser.add_argument("-g", "--genome_index", action="store", help="genome index", required=True)
parser.add_argument("-w", "--window_size", action="store", type= int, help="window size (bp)", required=True)
parser.add_argument("-o", "--output", action="store", help="output tab delimited file name", required=True)
args = parser.parse_args()

"""
# This script uses gwas output, a user defined window (bp), and gene annotations to profile genes around candidate causal SNPs
# input data
## input gwas output file
### format: SNP_name\tchromosome\tbasepair\tpvalue
## annotation file
### example format: gene_name\tgene_homolog_name\tgene_homolog_symbol\tgene_definition\tgene_homolog_tf_family\tgene_tf_family
## bed file
### format: chr\tstart\tend\tgene_name\tscore\tstrand
## genome index
### format: chr\tchromosome_length -> see samtools faidx
## output
### format: see annotation file format

"""

# read tab delimited file into pandas dataframe
snp_chr_bp_pvalue_df = pd.read_csv(args.input, sep = '\t')

# filter SNPs by desired pvalue 
snp_chr_bp_pvaluefilt_df = snp_chr_bp_pvalue_df[snp_chr_bp_pvalue_df["P"] < args.pvalue]

# adjust the start coordinate as bed files use a 0-based index
snp_chr_bp_pvaluefilt_df['start'] = snp_chr_bp_pvaluefilt_df['BP'] - 1

# reorder columns
snp_chr_bp_pvaluefilt_df_clean = snp_chr_bp_pvaluefilt_df[['CHR','start','BP', 'SNP']]

# change dtype of 'BP' and 'CHR' to string
snp_chr_bp_pvaluefilt_df_clean[['CHR','BP']] = snp_chr_bp_pvaluefilt_df_clean[['CHR','BP']].astype(str)

# concatenate columns with '-' delimeter
snp_chr_bp_pvaluefilt_df_clean['chr_snp_name'] = snp_chr_bp_pvaluefilt_df_clean[['CHR','BP','SNP']].apply(lambda x: '-'.join(x), axis=1)

# sort columns by chromosome and bp start
snp_chr_bp_pvaluefilt_df_clean.sort_values(['CHR', 'start'], inplace=True)

# reorder columns to bed format
snp_chr_bp_pvaluefilt_df_clean_fin = snp_chr_bp_pvaluefilt_df_clean[['CHR','start','BP', 'chr_snp_name']]

# convert pandas dataframe to string
snp_chr_bp_pvaluefilt_df_clean_fin_txt = pd.DataFrame(snp_chr_bp_pvaluefilt_df_clean_fin).to_csv(index = False, header = False, sep = '\t')

# convert to bedtools object
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool = pybedtools.BedTool(snp_chr_bp_pvaluefilt_df_clean_fin_txt, from_string = True)

# add upstream and downstream flanking regions
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_slop = snp_chr_bp_pvaluefilt_df_clean_fin_bedtool.slop(b=args.window_size, g=args.genome_index)

# intersect slop bed file with gene bed file
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_intersect = snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_slop.intersect(args.bed_file, wb = True, nonamecheck = True)

# read in intersect df
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_intersect_df = pd.read_csv(snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_intersect.fn, sep='\t', header=None)

# drop rows with duplicate genes, keeping the first genes
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates = snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_intersect_df.drop_duplicates(subset=[7], keep='first')

# read in tab delimited gene annotations as a pandas dataframe
annotation_df = pd.read_csv(args.annotation_file, sep='\t')

# rename gene column to the gene column name present in the annotation df
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates.rename(columns={7:list(annotation_df.columns.values)[0]}, inplace=True)

# merge gene dataframe and annotation dataframe
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates_merge = snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates.merge(annotation_df, how='left', on=[list(annotation_df.columns.values)[0]])

# subset columns
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates_merge_subset = snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates_merge[list(annotation_df.columns.values)]

# write out file
snp_chr_bp_pvaluefilt_df_clean_fin_bedtool_dropduplicates_merge_subset.to_csv(args.output, sep='\t', index=False)

