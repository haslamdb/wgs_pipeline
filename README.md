# WGS Pipeline

A comprehensive Python package for analyzing whole genome sequencing (WGS) data from bacterial samples. This package contains a structured Python package that handles file preparation, read trimming, assembly, quality assessment, annotation, and phylogenetic analysis.

## Features

### Core Features
- File preparation and organization
- Quality control with Trimmomatic and FastQC
- Genome assembly with Unicycler (Illumina-only or Hybrid)
- Assembly quality assessment with QUAST
- Coverage analysis with BBMap
- MLST typing

### Assembly & QC
- **Hybrid Assembly**: Unicycler with Illumina + Nanopore reads
- **Nanopore QC**: NanoPlot for long read quality assessment
- **Assembly QC**: QUAST for comprehensive assembly statistics

### Annotation
- **Prokka**: Traditional annotation (legacy)
- **Bakta**: Modern, comprehensive genome annotation

### Comparative Genomics & Phylogeny
- **Panaroo**: Pangenome analysis
- **Parsnp**: Core genome alignment and ML phylogenetic trees
- **Snippy**: Variant calling against reference genome

### Mutation Analysis
- **breseq**: Bacterial resequencing for mutation detection
  - SNPs, indels, structural variants
  - Mixed population detection
  - Supports Illumina, Nanopore, or hybrid consensus mode
- **gdtools**: Multi-sample mutation comparison

### Plasmid Analysis
- **MOB-suite**: Plasmid identification, reconstruction, and typing
  - Identifies plasmid contigs
  - Determines mobility type (conjugative/mobilizable)
  - Types incompatibility groups
  - Reconstructs plasmid sequences

### Analysis Outputs
- SNP matrices and distance matrices
- Phylogenetic trees (Newick format)
- HTML reports (breseq, QUAST)
- Comprehensive annotation files (GFF, GenBank, FASTA)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/wgs_pipeline.git
cd wgs_pipeline

# Install the package
pip install -e .
```

## Dependencies

This package requires several bioinformatics tools organized into separate conda environments for compatibility.

### Conda Environments Setup

**Status: ✅ All environments are installed and configured on this system!**

#### 1. Main Assembly Environment (unicycler-env) ✅
**Status:** INSTALLED
```bash
# Already installed! Contains:
# - unicycler (0.5.1)
# - trimmomatic (0.40)
# - fastqc (0.12.1)
# - bbmap (39.01)
```

To recreate on another system:
```bash
conda create -n unicycler-env -c bioconda -c conda-forge \
    unicycler trimmomatic fastqc bbmap \
    -y
```

**Tools:**
- Unicycler (hybrid/Illumina assembly)
- Trimmomatic (read trimming)
- FastQC (quality control)
- BBMap (coverage analysis)

#### 2. Nanopore QC Environment (nanopore-qc) ✅
**Status:** INSTALLED (separate environment)
```bash
# Already installed! Contains:
# - nanoplot (1.46.1)
# - nanostat (1.6.0)
```

**Tools:**
- NanoPlot (Nanopore read QC)
- NanoStat (Nanopore statistics)

**Note:** Use this environment specifically for `--step nanoplot`

#### 3. Assembly QC and Comparative Genomics (assembly-qc) ✅
**Status:** INSTALLED
```bash
# Already installed! Contains:
# - quast (5.2.0)
# - parsnp (1.5.6)
# - harvesttools (1.2)
# - fasttree (2.2.0)
```

To recreate on another system:
```bash
conda create -n assembly-qc -c bioconda -c conda-forge \
    quast parsnp harvesttools fasttree \
    -y
```

**Tools:**
- QUAST (assembly quality assessment)
- Parsnp (core genome alignment)
- HarvestTools (format conversion)
- FastTree (phylogenetic trees)

#### 4. Mutation Analysis (breseq) ✅
**Status:** INSTALLED
```bash
# Already installed! Contains:
# - breseq (0.39.0)
# - gdtools (included with breseq)
```

To recreate on another system:
```bash
conda create -n breseq -c bioconda breseq -y
```

**Tools:**
- breseq (bacterial resequencing)
- gdtools (mutation comparison)

#### 5. Modern Annotation (bakta) ✅
**Status:** INSTALLED
```bash
# Already installed! Contains:
# - bakta (1.11.4)
```

To recreate on another system:
```bash
conda create -n bakta -c bioconda -c conda-forge bakta -y
```

**Database Setup:**
```bash
# Download Bakta database (required, ~3GB for light version)
conda activate bakta
bakta_db download --output ~/bakta_db --type light

# Or full database (~30 GB)
bakta_db download --output ~/bakta_db --type full
```

**Tools:**
- Bakta (genome annotation)

#### 6. Plasmid Analysis (mobsuite) ✅
**Status:** INSTALLED
```bash
# Already installed! Contains:
# - mob_suite (3.1.9)
```

To recreate on another system:
```bash
conda create -n mobsuite -c bioconda -c conda-forge mob_suite -y
```

**Tools:**
- MOB-suite (plasmid identification and typing)

#### 7. Legacy Tools (wgs-legacy) ✅
**Status:** INSTALLED
```bash
# Already installed! Contains:
# - prokka (1.13)
# - mlst (2.11)
# - snippy (4.0.2)
# - bcftools (1.9)
# - fasttree (2.2.0)
```

To recreate on another system:
```bash
conda create -n wgs-legacy -c bioconda -c conda-forge \
    prokka mlst snippy bcftools fasttree \
    -y
```

**Tools:**
- MLST (sequence typing)
- Prokka (annotation - legacy)
- Snippy (variant calling)
- bcftools (VCF processing)
- FastTree (phylogenetic trees)

**Note:** Panaroo was excluded due to dependency conflicts. If needed, install separately.

### Quick Environment Check
```bash
# List all WGS-related environments
conda env list | grep -E "unicycler|assembly|breseq|bakta|mobsuite|wgs|nanopore"

# Expected output:
# assembly-qc
# bakta
# breseq
# mobsuite
# nanopore-qc
# unicycler-env
# wgs-legacy
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

### Workflow Examples

#### Example 1: Hybrid Assembly Workflow (Illumina + Nanopore)
```bash
# 1. Nanopore QC (use nanopore-qc environment)
conda activate nanopore-qc
wgs_pipeline --data_dir /path/to/data --ref_genome ref.fasta \
             --ref_gff ref.gff --ref_protein ref.faa \
             --sample_file samples.txt --step nanoplot

# 2. Hybrid Assembly (use unicycler-env)
conda activate unicycler-env
wgs_pipeline --data_dir /path/to/data --ref_genome ref.fasta \
             --ref_gff ref.gff --ref_protein ref.faa \
             --sample_file samples.txt --step hybrid_unicycler

# 3. Assembly QC (use assembly-qc environment)
conda activate assembly-qc
wgs_pipeline --data_dir /path/to/data --ref_genome ref.fasta \
             --ref_gff ref.gff --ref_protein ref.faa \
             --sample_file samples.txt --step quast

# 4. Core genome phylogeny (stay in assembly-qc)
wgs_pipeline --data_dir /path/to/data --ref_genome ref.fasta \
             --ref_gff ref.gff --ref_protein ref.faa \
             --sample_file samples.txt --step parsnp
```

#### Example 2: Mutation Detection Workflow
```bash
conda activate breseq

# Run breseq mutation analysis (Illumina mode)
wgs_pipeline --data_dir /path/to/data --ref_genome ref.fasta \
             --ref_gff ref.gff --ref_protein ref.faa \
             --sample_file samples.txt \
             --step breseq --breseq_mode illumina

# Compare mutations across samples
wgs_pipeline --data_dir /path/to/data --ref_genome ref.fasta \
             --ref_gff ref.gff --ref_protein ref.faa \
             --sample_file samples.txt --step gdtools
```

#### Example 3: Complete E. coli Analysis Pipeline
```bash
# Phase 1: Assembly (unicycler-env)
conda activate unicycler-env
wgs_pipeline --data_dir ~/Data/EcoliProject \
             --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff \
             --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt \
             --step hybrid_unicycler --threads 32

# Phase 2: QC & Phylogeny (assembly-qc)
conda activate assembly-qc
wgs_pipeline --data_dir ~/Data/EcoliProject --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt --step quast
wgs_pipeline --data_dir ~/Data/EcoliProject --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt --step parsnp

# Phase 3: Annotation (bakta)
conda activate bakta
wgs_pipeline --data_dir ~/Data/EcoliProject --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt --step bakta --bakta_db ~/bakta_db/db

# Phase 4: Plasmid Analysis (mobsuite)
conda activate mobsuite
wgs_pipeline --data_dir ~/Data/EcoliProject --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt --step mobsuite

# Phase 5: Mutation Analysis (breseq)
conda activate breseq
wgs_pipeline --data_dir ~/Data/EcoliProject --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt --step breseq --breseq_mode consensus
wgs_pipeline --data_dir ~/Data/EcoliProject --ref_genome ~/ref/ecoli_nissle.fasta \
             --ref_gff ~/ref/ecoli_nissle.gff --ref_protein ~/ref/ecoli_nissle.faa \
             --sample_file samples.txt --step gdtools
```

### Available Pipeline Steps

| Step | Conda Env | Description |
|------|-----------|-------------|
| `prepare` | unicycler-env | File preparation and merging |
| `trim` | unicycler-env | Trimmomatic read trimming |
| `fastqc` | unicycler-env | FastQC quality control |
| `nanoplot` | **nanopore-qc** | Nanopore read QC |
| `unicycler` | unicycler-env | Illumina-only assembly |
| `hybrid_unicycler` | unicycler-env | Hybrid assembly (Illumina + Nanopore) |
| `quast` | assembly-qc | Assembly quality assessment |
| `parsnp` | assembly-qc | Core genome alignment & phylogeny |
| `bbmap` | unicycler-env | Coverage analysis |
| `mlst` | wgs-legacy | MLST typing |
| `prokka` | wgs-legacy | Prokka annotation |
| `bakta` | bakta | Bakta annotation (requires --bakta_db) |
| `mobsuite` | mobsuite | Plasmid identification & typing |
| `panaroo` | wgs-legacy | Pangenome analysis (NOT INSTALLED) |
| `snippy` | wgs-legacy | Variant calling |
| `snippy_core` | wgs-legacy | Core genome SNP calling |
| `breseq` | breseq | Mutation detection (requires --breseq_mode) |
| `gdtools` | breseq | Multi-sample mutation comparison |
| `filter_snps` | wgs-legacy | SNP filtering |
| `distance_matrix` | wgs-legacy | SNP distance matrix |
| `clean` | any | Clean temporary files |

**Note:** Panaroo is not currently installed due to dependency conflicts. To use pangenome analysis, install Panaroo manually or skip this step.

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
