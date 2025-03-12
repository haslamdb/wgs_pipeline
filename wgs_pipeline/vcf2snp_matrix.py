#!/usr/bin/env python3
"""
vcf2SNPmatrix.py

This script parses a VCF file to create a SNP matrix and pairwise distance matrix.
It's intended to be used with the output from snippy-core.

Usage:
    python vcf2SNPmatrix.py -i core.vcf -o snp_matrix.csv

Example with snippy-core output:
    python vcf2SNPmatrix.py -i snippy_core/core.vcf -o snippy_core/snp_matrix.csv
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Convert VCF to SNP matrix and distance matrix')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', default='snp_matrix.csv', help='Output matrix CSV file')
    parser.add_argument('-d', '--distance', action='store_true', help='Calculate distance matrix')
    parser.add_argument('-f', '--filtered', action='store_true', help='Filter for biallelic SNPs only')
    parser.add_argument('-m', '--min_alt_freq', type=float, default=0.05, 
                        help='Minimum alternate allele frequency (default: 0.05)')
    
    return parser.parse_args()

def read_vcf_header(vcf_file):
    """Read the header of a VCF file to get sample names."""
    samples = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # Parse the header line to get sample names
                fields = line.strip().split('\t')
                samples = fields[9:]  # Sample names start at column 10
                break
    
    return samples

def parse_vcf(vcf_file, filter_biallelic=True, min_alt_freq=0.05):
    """
    Parse VCF file and extract the SNP data.
    
    Args:
        vcf_file: Path to the VCF file
        filter_biallelic: Whether to filter for biallelic SNPs only
        min_alt_freq: Minimum alternate allele frequency
        
    Returns:
        positions: List of SNP positions (CHROM_POS)
        ref_alleles: List of reference alleles
        alt_alleles: List of alternate alleles
        genotypes: Dictionary mapping sample names to lists of genotypes
    """
    # Get sample names from the header
    samples = read_vcf_header(vcf_file)
    
    positions = []
    ref_alleles = []
    alt_alleles = []
    genotypes = {sample: [] for sample in samples}
    
    # Parse the VCF file
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            # Filter for SNPs only (exclude indels)
            if len(ref) > 1 or any(len(a) > 1 for a in alt.split(',')):
                continue
            
            # Filter for biallelic SNPs if requested
            if filter_biallelic and ',' in alt:
                continue
            
            # Calculate alternate allele frequency
            sample_data = fields[9:]
            genotype_counts = defaultdict(int)
            
            for data in sample_data:
                gt = data.split(':')[0]  # Extract genotype field (GT)
                genotype_counts[gt] += 1
            
            # Skip sites where the alternate allele frequency is too low
            if len(samples) > 0:
                alt_count = sum(genotype_counts[gt] for gt in genotype_counts if gt != '0/0' and gt != '.')
                alt_freq = alt_count / len(samples)
                
                if alt_freq < min_alt_freq:
                    continue
            
            # Add position and alleles
            positions.append(f"{chrom}_{pos}")
            ref_alleles.append(ref)
            alt_alleles.append(alt)
            
            # Add genotypes for each sample
            for i, sample in enumerate(samples):
                if i < len(sample_data):
                    gt = sample_data[i].split(':')[0]
                    genotypes[sample].append(gt)
                else:
                    genotypes[sample].append('.')
    
    return positions, ref_alleles, alt_alleles, genotypes

def create_snp_matrix(positions, ref_alleles, alt_alleles, genotypes):
    """
    Create a SNP matrix from the parsed VCF data.
    
    Args:
        positions: List of SNP positions (CHROM_POS)
        ref_alleles: List of reference alleles
        alt_alleles: List of alternate alleles
        genotypes: Dictionary mapping sample names to lists of genotypes
        
    Returns:
        pandas.DataFrame: SNP matrix
    """
    # Create a DataFrame with positions as index
    matrix_data = {
        'REF': ref_alleles,
        'ALT': alt_alleles
    }
    
    # Add sample genotypes
    for sample, gts in genotypes.items():
        matrix_data[sample] = gts
    
    # Create DataFrame
    snp_matrix = pd.DataFrame(matrix_data, index=positions)
    
    return snp_matrix

def create_binary_matrix(snp_matrix):
    """
    Convert the SNP matrix to a binary matrix (0 for reference, 1 for alternate).
    
    Args:
        snp_matrix: SNP matrix from create_snp_matrix
        
    Returns:
        pandas.DataFrame: Binary matrix
    """
    # Create a copy of the original matrix
    binary_matrix = snp_matrix.copy()
    
    # Get sample columns (exclude REF and ALT)
    sample_cols = [col for col in binary_matrix.columns if col not in ['REF', 'ALT']]
    
    # Convert genotypes to binary values
    for col in sample_cols:
        binary_matrix[col] = binary_matrix[col].apply(
            lambda gt: 0 if gt == '0/0' or gt == '.' else 1
        )
    
    return binary_matrix

def create_distance_matrix(binary_matrix):
    """
    Create a pairwise distance matrix from the binary matrix.
    
    Args:
        binary_matrix: Binary matrix from create_binary_matrix
        
    Returns:
        pandas.DataFrame: Distance matrix
    """
    # Get sample columns (exclude REF and ALT)
    sample_cols = [col for col in binary_matrix.columns if col not in ['REF', 'ALT']]
    
    # Create an empty distance matrix
    distance_matrix = pd.DataFrame(index=sample_cols, columns=sample_cols)
    
    # Fill the distance matrix with pairwise SNP distances
    for i, sample1 in enumerate(sample_cols):
        for j, sample2 in enumerate(sample_cols):
            if i == j:
                distance_matrix.loc[sample1, sample2] = 0
            else:
                # Count SNP differences
                diff_count = (binary_matrix[sample1] != binary_matrix[sample2]).sum()
                distance_matrix.loc[sample1, sample2] = diff_count
    
    return distance_matrix

def main():
    """Main function."""
    args = parse_args()
    
    # Check if the input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        sys.exit(1)
    
    print(f"Parsing VCF file: {args.input}")
    positions, ref_alleles, alt_alleles, genotypes = parse_vcf(
        args.input, 
        filter_biallelic=args.filtered,
        min_alt_freq=args.min_alt_freq
    )
    
    print(f"Found {len(positions)} SNPs")
    
    # Create SNP matrix
    snp_matrix = create_snp_matrix(positions, ref_alleles, alt_alleles, genotypes)
    
    # Create binary matrix
    binary_matrix = create_binary_matrix(snp_matrix)
    
    # Save SNP matrix to CSV
    output_file = args.output
    binary_matrix.to_csv(output_file)
    print(f"Saved SNP matrix to: {output_file}")
    
    # Create and save distance matrix if requested
    if args.distance:
        distance_matrix = create_distance_matrix(binary_matrix)
        
        # Save distance matrix to CSV
        distance_output = output_file.replace('.csv', '_distance.csv')
        distance_matrix.to_csv(distance_output)
        print(f"Saved distance matrix to: {distance_output}")
    
    print("Done!")

if __name__ == "__main__":
    main()
