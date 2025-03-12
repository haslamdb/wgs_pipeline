#!/usr/bin/env python3
"""
WGS Analysis Pipeline for Bacterial Genomes

This package converts the bash script workflow for analyzing whole genome sequencing data
from bacterial samples into a structured Python package. It handles file preparation,
read trimming, assembly, quality assessment, annotation, and phylogenetic analysis.

Example usage:
    from wgs_pipeline import WGSPipeline
    
    pipeline = WGSPipeline(
        data_dir="~/Data/WGSData/AndreaData",
        reference_genome="/home/david/Databases/BacterialDatabases/GCF_000009205.2_ASM920v2_genomic.fasta",
        reference_gff="/home/david/Databases/BacterialDatabases/GCF_000009205.2_ASM920v2_genomic.genomic.gff",
        reference_protein="/home/david/Databases/BacterialDatabases/GGCF_000210435.1_ASM21043v1_protein.faa"
    )
    
    # Run the full pipeline
    pipeline.run_pipeline("AndreaMupRFilesRevised.txt")
    
    # Or run specific steps
    pipeline.prepare_files("AndreaMupRFilesRevised.txt")
    pipeline.run_trimmomatic("AndreaMupRFilesRevised.txt")
    # ...etc
"""

import os
import sys
import subprocess
import glob
import pandas as pd
import allel
import shutil
from pathlib import Path
import logging
from typing import List, Dict, Optional, Union, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("wgs_pipeline.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("wgs_pipeline")

class WGSPipeline:
    """
    A class to run a whole genome sequencing analysis pipeline for bacterial genomes.
    
    This class implements methods that correspond to steps in the original bash script,
    including file preparation, read trimming, assembly, and more advanced analyses.
    """
    
    def __init__(self, 
                 data_dir: str,
                 reference_genome: str,
                 reference_gff: str,
                 reference_protein: str,
                 threads: int = 48):
        """
        Initialize the WGS pipeline with paths to important directories and reference files.
        
        Args:
            data_dir: Base directory containing sequencing data
            reference_genome: Path to reference genome FASTA
            reference_gff: Path to reference genome GFF annotation
            reference_protein: Path to reference protein FASTA
            threads: Number of threads to use for parallel processing
        """
        # Expand paths
        self.data_dir = os.path.expanduser(data_dir)
        self.reference_genome = os.path.expanduser(reference_genome)
        self.reference_gff = os.path.expanduser(reference_gff)
        self.reference_protein = os.path.expanduser(reference_protein)
        
        # Create directories if they don't exist
        self.gff_dir = os.path.join(os.path.dirname(self.data_dir), "GFFFiles")
        self.snippy_dir = os.path.join(os.path.dirname(self.data_dir), "SnippyVCF")
        
        os.makedirs(self.gff_dir, exist_ok=True)
        os.makedirs(self.snippy_dir, exist_ok=True)
        
        # Set threads for parallel processing
        self.threads = threads
        
        logger.info(f"Initialized WGS pipeline with data directory: {self.data_dir}")
        logger.info(f"Reference genome: {self.reference_genome}")
        
    def _run_command(self, command: str, cwd: Optional[str] = None) -> Tuple[int, str, str]:
        """
        Run a shell command and log the output.
        
        Args:
            command: Command to run
            cwd: Directory to run the command in (optional)
            
        Returns:
            Tuple containing (return_code, stdout, stderr)
        """
        logger.info(f"Running command: {command}")
        
        process = subprocess.Popen(
            command, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            cwd=cwd,
            text=True
        )
        
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Command failed with code {process.returncode}")
            logger.error(f"stderr: {stderr}")
        else:
            logger.info(f"Command completed successfully")
            
        return process.returncode, stdout, stderr
    
    def _read_sample_list(self, sample_file: str) -> List[str]:
        """
        Read the sample list file and return a list of sample IDs.
        
        Args:
            sample_file: Path to file containing sample names
            
        Returns:
            List of sample IDs
        """
        file_path = os.path.join(self.data_dir, sample_file)
        
        with open(file_path, 'r') as f:
            samples = [line.strip() for line in f if line.strip()]
            
        logger.info(f"Read {len(samples)} samples from {file_path}")
        return samples
    
    def prepare_files(self, sample_file: str, lane8_dir: str = "Lane8") -> None:
        """
        Prepare files by merging lanes and renaming files.
        
        Args:
            sample_file: File containing sample names
            lane8_dir: Directory containing Lane 8 reads
        """
        samples = self._read_sample_list(sample_file)
        lane8_path = os.path.join(self.data_dir, lane8_dir)
        
        logger.info(f"Preparing files for {len(samples)} samples")
        logger.info(f"Lane 8 directory: {lane8_path}")
        
        for sample in samples:
            # Merge R1 files
            r1_lane8 = os.path.join(lane8_path, f"{sample}_R1.fastq.gz")
            r1_output = os.path.join(self.data_dir, f"{sample}_R1.fastq.gz")
            
            if os.path.exists(r1_lane8) and os.path.exists(r1_output):
                self._run_command(f"cat {r1_lane8} >> {r1_output}")
                
            # Merge R2 files
            r2_lane8 = os.path.join(lane8_path, f"{sample}_R2.fastq.gz")
            r2_output = os.path.join(self.data_dir, f"{sample}_R2.fastq.gz")
            
            if os.path.exists(r2_lane8) and os.path.exists(r2_output):
                self._run_command(f"cat {r2_lane8} >> {r2_output}")
                
        logger.info("File preparation completed")
    
    def run_trimmomatic(self, sample_file: str) -> None:
        """
        Run Trimmomatic on all samples to trim low-quality reads.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        
        logger.info(f"Running Trimmomatic on {len(samples)} samples")
        
        for sample in samples:
            input_r1 = os.path.join(self.data_dir, f"{sample}_R1.fastq.gz")
            input_r2 = os.path.join(self.data_dir, f"{sample}_R2.fastq.gz")
            
            output_r1 = os.path.join(self.data_dir, f"trimmed_{sample}_R1.fastq.gz")
            output_r2 = os.path.join(self.data_dir, f"trimmed_{sample}_R2.fastq.gz")
            
            unpaired_r1 = os.path.join(self.data_dir, f"unpaired_{sample}_R1.fastq.gz")
            unpaired_r2 = os.path.join(self.data_dir, f"unpaired_{sample}_R2.fastq.gz")
            
            cmd = (
                f"java -jar /usr/local/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE "
                f"{input_r1} {input_r2} "
                f"{output_r1} {unpaired_r1} "
                f"{output_r2} {unpaired_r2} "
                f"LEADING:5 TRAILING:5 MINLEN:60"
            )
            
            self._run_command(cmd)
            
            # Remove unpaired reads
            self._run_command(f"rm {unpaired_r1} {unpaired_r2}")
            
        logger.info("Trimmomatic completed for all samples")
    
    def run_fastqc(self, sample_file: str) -> None:
        """
        Run FastQC on trimmed reads.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        
        logger.info(f"Running FastQC on {len(samples)} samples")
        
        for sample in samples:
            trimmed_r1 = os.path.join(self.data_dir, f"trimmed_{sample}_R1.fastq.gz")
            trimmed_r2 = os.path.join(self.data_dir, f"trimmed_{sample}_R2.fastq.gz")
            
            # Create sample directory if it doesn't exist
            sample_dir = os.path.join(self.data_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)
            
            # Run FastQC
            cmd = f"fastqc {trimmed_r1} {trimmed_r2} --extract"
            self._run_command(cmd)
            
            # Move FastQC output to sample directory
            self._run_command(f"mv trimmed_{sample}_R1_fastqc/ {sample_dir}")
            self._run_command(f"mv trimmed_{sample}_R2_fastqc/ {sample_dir}")
            self._run_command(f"mv trimmed_{sample}_R1_fastqc.html {sample_dir}")
            self._run_command(f"mv trimmed_{sample}_R2_fastqc.html {sample_dir}")
            
        logger.info("FastQC completed for all samples")
    
    def run_unicycler(self, sample_file: str) -> None:
        """
        Assemble reads into contigs using Unicycler.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        
        logger.info(f"Running Unicycler assembly on {len(samples)} samples")
        
        for sample in samples:
            trimmed_r1 = os.path.join(self.data_dir, f"trimmed_{sample}_R1.fastq.gz")
            trimmed_r2 = os.path.join(self.data_dir, f"trimmed_{sample}_R2.fastq.gz")
            
            # Create sample directory if it doesn't exist
            sample_dir = os.path.join(self.data_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)
            
            # Create output directory for unicycler
            unicycler_dir = os.path.join(sample_dir, f"{sample}_Unicycler")
            
            # Run Unicycler
            cmd = f"unicycler -1 {trimmed_r1} -2 {trimmed_r2} -t {self.threads} -o {unicycler_dir}"
            self._run_command(cmd)
            
            # Copy assembly to the main WGS data directory
            wgs_data_dir = os.path.dirname(self.data_dir)
            sample_wgs_dir = os.path.join(wgs_data_dir, sample)
            os.makedirs(sample_wgs_dir, exist_ok=True)
            
            self._run_command(f"cp {unicycler_dir}/assembly.fasta {sample_wgs_dir}/{sample}.fasta")
            self._run_command(f"cp {unicycler_dir}/assembly.gfa {sample_wgs_dir}/{sample}.gfa")
            
            # Also copy to WGS root directory for convenience
            self._run_command(f"cp {unicycler_dir}/assembly.fasta {wgs_data_dir}/{sample}.fasta")
            self._run_command(f"cp {unicycler_dir}/assembly.gfa {wgs_data_dir}/{sample}.gfa")
            
        logger.info("Unicycler assembly completed for all samples")
    
    def run_quast(self, sample_file: str) -> None:
        """
        Run Quast to get assembly statistics.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        wgs_data_dir = os.path.dirname(self.data_dir)
        
        logger.info(f"Running Quast on {len(samples)} samples")
        
        for sample in samples:
            # Construct the file name
            fasta_file = os.path.join(wgs_data_dir, f"{sample}.fasta")
            
            if not os.path.exists(fasta_file):
                logger.warning(f"Fasta file not found: {fasta_file}, skipping Quast")
                continue
            
            # Run QUAST
            cmd = (
                f"quast.py {fasta_file} "
                f"-o {sample}_Quast "
                f"-r {self.reference_genome} "
                f"-g {self.reference_gff}"
            )
            
            self._run_command(cmd, cwd=wgs_data_dir)
            
            # Copy the report
            self._run_command(
                f"cp {sample}_Quast/report.tsv {sample}_Quast_report.tsv",
                cwd=wgs_data_dir
            )
            
        logger.info("Quast completed for all samples")
    
    def run_bbmap(self, sample_file: str) -> None:
        """
        Run BBMap to get coverage statistics.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        
        logger.info(f"Running BBMap on {len(samples)} samples")
        
        for sample in samples:
            trimmed_r1 = os.path.join(self.data_dir, f"trimmed_{sample}_R1.fastq.gz")
            trimmed_r2 = os.path.join(self.data_dir, f"trimmed_{sample}_R2.fastq.gz")
            
            # Create sample directory if it doesn't exist
            sample_dir = os.path.join(self.data_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)
            
            # Run BBMap
            cmd = (
                f"bbmap.sh ref={self.reference_genome} "
                f"in={trimmed_r1} in2={trimmed_r2} "
                f"out={sample}.sam outu={sample}_unmapped.fastq.gz "
                f"covstats={sample}_constats.txt bincov={sample}_bincov.txt"
            )
            
            self._run_command(cmd, cwd=self.data_dir)
            
            # Move output files to sample directory
            self._run_command(f"mv {sample}.sam {sample_dir}/{sample}_bbmap.sam", cwd=self.data_dir)
            self._run_command(f"mv {sample}_unmapped.fastq.gz {sample_dir}/{sample}_unmapped.fastq.gz", cwd=self.data_dir)
            self._run_command(f"mv {sample}_constats.txt {sample_dir}/{sample}_constats.txt", cwd=self.data_dir)
            self._run_command(f"mv {sample}_bincov.txt {sample_dir}/{sample}_bincov.txt", cwd=self.data_dir)
            
        logger.info("BBMap completed for all samples")
    
    def run_mlst(self, sample_file: str) -> None:
        """
        Run MLST to determine sequence type.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        wgs_data_dir = os.path.dirname(self.data_dir)
        
        logger.info(f"Running MLST on {len(samples)} samples")
        
        # Create output file
        mlst_output = os.path.join(wgs_data_dir, "MLST_Table.tsv")
        
        for sample in samples:
            fasta_file = os.path.join(wgs_data_dir, f"{sample}.fasta")
            
            if not os.path.exists(fasta_file):
                logger.warning(f"Fasta file not found: {fasta_file}, skipping MLST")
                continue
            
            # Run MLST
            cmd = f"mlst {fasta_file} >> {mlst_output}"
            self._run_command(cmd, cwd=wgs_data_dir)
            
        logger.info("MLST completed for all samples")
    
    def run_prokka(self, sample_file: str) -> None:
        """
        Run Prokka to annotate genomes.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        wgs_data_dir = os.path.dirname(self.data_dir)
        
        logger.info(f"Running Prokka on {len(samples)} samples")
        
        for sample in samples:
            fasta_file = os.path.join(wgs_data_dir, f"{sample}.fasta")
            
            if not os.path.exists(fasta_file):
                logger.warning(f"Fasta file not found: {fasta_file}, skipping Prokka")
                continue
            
            # Create output directory
            sample_dir = os.path.join(wgs_data_dir, sample)
            prokka_dir = os.path.join(sample_dir, f"{sample}_Prokka")
            os.makedirs(sample_dir, exist_ok=True)
            
            # Run Prokka
            cmd = (
                f"prokka -mincontiglen 500 --force --cpus {self.threads} --rnammer "
                f"--outdir {prokka_dir} --prefix {sample} --proteins {self.reference_protein} {fasta_file}"
            )
            
            self._run_command(cmd, cwd=wgs_data_dir)
            
            # Move GFF file to GFF directory for Panaroo
            self._run_command(f"cp {prokka_dir}/*.gff {self.gff_dir}/{sample}.gff", cwd=wgs_data_dir)
            
        logger.info("Prokka completed for all samples")
    
    def run_panaroo(self) -> None:
        """
        Run Panaroo to generate a pangenome analysis.
        """
        logger.info("Running Panaroo for pangenome analysis")
        
        # Get all GFF files
        gff_files = glob.glob(os.path.join(self.gff_dir, "*.gff"))
        
        if not gff_files:
            logger.warning("No GFF files found in GFF directory, skipping Panaroo")
            return
        
        # Run Panaroo
        output_dir = os.path.join(self.gff_dir, "Panaroo_Out")
        cmd = f"panaroo -i {self.gff_dir}/*.gff -o {output_dir}"
        
        self._run_command(cmd)
        
        # Run FastTree on the output
        fasttree_cmd = (
            f"FastTree -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -nosupport -intree "
            f"{output_dir}/combined_DNA_CDS.fasta > {output_dir}/CombinedAndreaMupR.tre"
        )
        
        self._run_command(fasttree_cmd)
        
        logger.info("Panaroo and FastTree completed")
    
    def run_snippy(self, sample_file: str) -> None:
        """
        Run Snippy for variant calling.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        wgs_data_dir = os.path.dirname(self.data_dir)
        
        logger.info(f"Running Snippy on {len(samples)} samples")
        
        for sample in samples:
            trimmed_r1 = os.path.join(self.data_dir, f"trimmed_{sample}_R1.fastq.gz")
            trimmed_r2 = os.path.join(self.data_dir, f"trimmed_{sample}_R2.fastq.gz")
            
            # Create sample directory if it doesn't exist
            sample_dir = os.path.join(wgs_data_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)
            
            # Create output directory for Snippy
            snippy_dir = os.path.join(sample_dir, f"{sample}_Snippy")
            
            # Run Snippy
            cmd = (
                f"snippy --cpus {self.threads} --outdir {snippy_dir} "
                f"--ref {self.reference_genome} "
                f"--R1 {trimmed_r1} --R2 {trimmed_r2}"
            )
            
            self._run_command(cmd, cwd=wgs_data_dir)
            
            # Copy output to SnippyVCF directory
            self._run_command(f"cp -R {snippy_dir} {self.snippy_dir}")
            
        logger.info("Snippy completed for all samples")
    
    def run_snippy_core(self) -> None:
        """
        Run snippy-core to create a core genome alignment.
        """
        logger.info("Running snippy-core for core genome alignment")
        
        # Get all Snippy directories
        snippy_dirs = glob.glob(os.path.join(self.snippy_dir, "*_Snippy"))
        
        if not snippy_dirs:
            logger.warning("No Snippy directories found, skipping snippy-core")
            return
        
        # Run snippy-core
        cmd = f"snippy-core --ref {self.reference_genome} {self.snippy_dir}/*_Snippy"
        self._run_command(cmd, cwd=self.snippy_dir)
        
        # Run FastTree on the core alignment
        fasttree_cmd = (
            f"FastTree -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -nosupport "
            f"-intree {self.snippy_dir}/core.full.aln > {self.snippy_dir}/CoreAndreaMupR.tre"
        )
        self._run_command(fasttree_cmd)
        
        logger.info("snippy-core and FastTree completed")
    
    def filter_snps(self) -> None:
        """
        Filter SNPs from the core.vcf file and create a SNP matrix.
        """
        logger.info("Filtering SNPs and creating SNP matrix")
        
        # Check if core.vcf exists
        core_vcf = os.path.join(self.snippy_dir, "core.vcf")
        
        if not os.path.exists(core_vcf):
            logger.warning("core.vcf not found, skipping SNP filtering")
            return
        
        # Filter SNPs using bcftools
        filtered_vcf = os.path.join(self.snippy_dir, "filtered_snps.vcf")
        cmd = f"bcftools view -v snps -m2 -M2 -q 0.05:minor -Oz -o {filtered_vcf} {core_vcf}"
        self._run_command(cmd)
        
        # Create SNP matrix
        snp_matrix = os.path.join(self.snippy_dir, "snp_matrix.txt")
        cmd = f"bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' {filtered_vcf} > {snp_matrix}"
        self._run_command(cmd)
        
        # Get sample names
        sample_names = os.path.join(self.snippy_dir, "samplenames.txt")
        cmd = f"bcftools query -l {core_vcf} > {sample_names}"
        self._run_command(cmd)
        
        logger.info("SNP filtering and matrix creation completed")
    
    def create_snp_distance_matrix(self) -> None:
        """
        Create a SNP distance matrix from the SNP matrix.
        """
        logger.info("Creating SNP distance matrix")
        
        # Check if SNP matrix exists
        snp_matrix_file = os.path.join(self.snippy_dir, "snp_matrix.txt")
        sample_names_file = os.path.join(self.snippy_dir, "samplenames.txt")
        
        if not os.path.exists(snp_matrix_file) or not os.path.exists(sample_names_file):
            logger.warning("SNP matrix or sample names file not found, skipping distance matrix creation")
            return
        
        try:
            # Read sample names
            with open(sample_names_file, 'r') as f:
                sample_names = [line.strip() for line in f if line.strip()]
            
            # Read SNP matrix
            snp_data = pd.read_csv(snp_matrix_file, sep='\t', header=None)
            
            # Convert genotypes to binary (presence/absence of alternate allele)
            genotype_data = snp_data.iloc[:, 4:]
            
            # Convert genotype strings to numeric values
            # '0/0' -> 0, '1/1' -> 1, etc.
            binary_matrix = genotype_data.applymap(lambda x: 0 if x == '0/0' else 1)
            
            # Set column names to sample names
            binary_matrix.columns = sample_names
            
            # Calculate pairwise distances
            num_samples = len(sample_names)
            distances = pd.DataFrame(index=sample_names, columns=sample_names)
            
            for i in range(num_samples):
                for j in range(num_samples):
                    # Calculate Hamming distance (number of SNP differences)
                    if i == j:
                        distances.iloc[i, j] = 0
                    else:
                        # Count differences
                        diff_count = (binary_matrix.iloc[:, i] != binary_matrix.iloc[:, j]).sum()
                        distances.iloc[i, j] = diff_count
            
            # Save distance matrix
            distances.to_csv(os.path.join(self.snippy_dir, "snp_distance_matrix.csv"))
            
            logger.info("SNP distance matrix created successfully")
            
        except Exception as e:
            logger.error(f"Error creating SNP distance matrix: {str(e)}")
    
    def clean_temp_files(self, sample_file: str) -> None:
        """
        Clean temporary files to save disk space.
        
        Args:
            sample_file: File containing sample names
        """
        samples = self._read_sample_list(sample_file)
        wgs_data_dir = os.path.dirname(self.data_dir)
        
        logger.info(f"Cleaning temporary files for {len(samples)} samples")
        
        for sample in samples:
            sample_dir = os.path.join(wgs_data_dir, sample)
            
            if not os.path.exists(sample_dir):
                continue
            
            # Remove large temporary files
            temp_files = [
                f"{sample}_bbmap.sam",
                "*spades_graph*",
                "*unmapped*",
                "*overlaps_removed*",
                "*bridges*",
                "*final_clean*"
            ]
            
            for pattern in temp_files:
                for file in glob.glob(os.path.join(sample_dir, pattern)):
                    os.remove(file)
            
            # Remove Snippy BAM files
            snippy_dir = os.path.join(sample_dir, f"{sample}_Snippy")
            if os.path.exists(snippy_dir):
                for file in glob.glob(os.path.join(snippy_dir, "*snps.bam")):
                    os.remove(file)
        
        logger.info("Temporary files cleaned")
    
    def run_pipeline(self, sample_file: str) -> None:
        """
        Run the complete pipeline from start to finish.
        
        Args:
            sample_file: File containing sample names
        """
        logger.info("Starting complete WGS pipeline")
        
        # Run all steps in sequence
        self.prepare_files(sample_file)
        self.run_trimmomatic(sample_file)
        self.run_fastqc(sample_file)
        self.run_unicycler(sample_file)
        self.run_quast(sample_file)
        self.run_bbmap(sample_file)
        self.run_mlst(sample_file)
        self.run_prokka(sample_file)
        self.run_panaroo()
        self.run_snippy(sample_file)
        self.run_snippy_core()
        self.filter_snps()
        self.create_snp_distance_matrix()
        self.clean_temp_files(sample_file)
        
        logger.info("Complete WGS pipeline finished")


# Command-line interface
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="WGS Analysis Pipeline")
    parser.add_argument("--data_dir", required=True, help="Directory containing sequencing data")
    parser.add_argument("--ref_genome", required=True, help="Path to reference genome FASTA")
    parser.add_argument("--ref_gff", required=True, help="Path to reference genome GFF annotation")
    parser.add_argument("--ref_protein", required=True, help="Path to reference protein FASTA")
    parser.add_argument("--sample_file", required=True, help="File containing sample names")
    parser.add_argument("--threads", type=int, default=48, help="Number of threads to use")
    parser.add_argument("--step", choices=[
        "all", "prepare", "trim", "fastqc", "unicycler", "quast", 
        "bbmap", "mlst", "prokka", "panaroo", "snippy", "snippy_core", "filter_snps", "distance_matrix", "clean"
    ], default="all", help="Pipeline step to run")
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = WGSPipeline(
        data_dir=args.data_dir,
        reference_genome=args.ref_genome,
        reference_gff=args.ref_gff,
        reference_protein=args.ref_protein,
        threads=args.threads
    )
    
    # Run selected pipeline step
    if args.step == "all":
        pipeline.run_pipeline(args.sample_file)
    elif args.step == "prepare":
        pipeline.prepare_files(args.sample_file)
    elif args.step == "trim":
        pipeline.run_trimmomatic(args.sample_file)
    elif args.step == "fastqc":
        pipeline.run_fastqc(args.sample_file)
    elif args.step == "unicycler":
        pipeline.run_unicycler(args.sample_file)
    elif args.step == "quast":
        pipeline.run_quast(args.sample_file)
    elif args.step == "bbmap":
        pipeline.run_bbmap(args.sample_file)
    elif args.step == "mlst":
        pipeline.run_mlst(args.sample_file)
    elif args.step == "prokka":
        pipeline.run_prokka(args.sample_file)
    elif args.step == "panaroo":
        pipeline.run_panaroo()
    elif args.step == "snippy":
        pipeline.run_snippy(args.sample_file)
    elif args.step == "snippy_core":
        pipeline.run_snippy_core()
    elif args.step == "filter_snps":
        pipeline.filter_snps()
    elif args.step == "distance_matrix":
        pipeline.create_snp_distance_matrix()
    elif args.step == "clean":
        pipeline.clean_temp_files(args.sample_file)