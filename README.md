# WGS Pipeline

A comprehensive Python package for analyzing whole genome sequencing (WGS) data from bacterial samples. This package contains a structured Python package that handles file preparation, read trimming, assembly, quality assessment, annotation, and phylogenetic analysis.

## Features

- File preparation and organization
- Quality control with Trimmomatic and FastQC
- Genome assembly with Unicycler
- Assembly quality assessment with Quast
- Coverage analysis with BBMap
- MLST typing
- Genome annotation with Prokka
- Pangenome analysis with Panaroo
- Variant calling with Snippy
- SNP analysis and distance matrix generation

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/wgs_pipeline.git
cd wgs_pipeline

# Install the package
pip install -e .
```

## Dependencies

This package requires several bioinformatics tools to be installed and available in your PATH:

- Trimmomatic
- FastQC
- Unicycler
- Quast
- BBMap
- CheckM
- MUMmer (dnadiff)
- MLST
- Prokka
- Panaroo
- Snippy
- FastTree
- bcftools

Most of these can be installed via conda:

```bash
conda create -n wgs_env
conda activate wgs_env
conda install -c bioconda trimmomatic fastqc unicycler quast bbmap mlst prokka panaroo snippy fasttree bcftools
```

## Usage

### Command Line Interface

```bash
# Run the full pipeline
wgs_pipeline --data_dir /path/to/raw/data \
             --ref_genome /path/to/reference.fasta \
             --ref_gff /path/to/reference.gff \
             --ref_protein /path/to/reference_proteins.faa \
             --sample_file /path/to/filelist.txt

# Run a specific step
wgs_pipeline --data_dir /path/to/raw/data \
             --ref_genome /path/to/reference.fasta \
             --ref_gff /path/to/reference.gff \
             --ref_protein /path/to/reference_proteins.faa \
             --sample_file /path/to/filelist.txt \
             --step unicycler
```

### Python API

```python
from wgs_pipeline import WGSPipeline

# Initialize the pipeline
pipeline = WGSPipeline(
    data_dir="/path/to/raw/data",
    reference_genome="/path/to/reference.fasta",
    reference_gff="/path/to/reference.gff",
    reference_protein="/path/to/reference_proteins.faa"
)

# Run the full pipeline
pipeline.run_pipeline("/path/to/filelist.txt")

# Or run specific steps
pipeline.prepare_files("/path/to/filelist.txt")
pipeline.run_trimmomatic("/path/to/filelist.txt")
pipeline.run_unicycler("/path/to/filelist.txt")
# ...etc
```

## VCF to SNP Matrix Utility

The package also includes a utility script for converting VCF files to SNP matrices:

```bash
vcf2snpmatrix -i snippy_core/core.vcf -o snippy_core/snp_matrix.csv -d
```

Options:
- `-i, --input`: Input VCF file (required)
- `-o, --output`: Output matrix CSV file (default: snp_matrix.csv)
- `-d, --distance`: Calculate distance matrix
- `-f, --filtered`: Filter for biallelic SNPs only
- `-m, --min_alt_freq`: Minimum alternate allele frequency (default: 0.05)

## Pipeline Structure

The pipeline follows these main steps:

1. **File Preparation**
   - Merge lanes and rename files

2. **Quality Control**
   - Trimmomatic: Trim low-quality reads
   - FastQC: Quality assessment of reads

3. **Assembly**
   - Unicycler: Assemble reads into contigs
   - Quast: Assembly quality assessment

4. **Coverage Analysis**
   - BBMap: Map reads to reference genome and calculate coverage

5. **Strain Typing**
   - MLST: Multi-Locus Sequence Typing

6. **Annotation**
   - Prokka: Genome annotation

7. **Pangenome Analysis**
   - Panaroo: Generate pangenomes and phylogenetic trees

8. **Variant Calling**
   - Snippy: Call variants against reference genome
   - snippy-core: Generate core genome alignment

9. **SNP Analysis**
   - Filter SNPs
   - Create SNP matrix and distance matrix

## File Organization

The pipeline organizes files in the following structure:

```
Data/
└── WGSData/
    ├── WGSData/                  # Raw data directory
    │   ├── trimmed_*.fastq.gz       # Trimmed reads
    │   └── [sample]/                # Sample directories
    │       ├── [sample]_Unicycler/  # Unicycler output
    │       ├── [sample]_Prokka/     # Prokka output
    │       └── [sample]_Snippy/     # Snippy output
    ├── GFFFiles/                    # GFF files for Panaroo
    │   └── Panaroo_Out/             # Panaroo output
    └── SnippyVCF/                   # Snippy output for core genome analysis
        ├── core.vcf                 # Core genome VCF
        ├── core.full.aln           # Core genome alignment
        └── snp_matrix.csv          # SNP matrix
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
