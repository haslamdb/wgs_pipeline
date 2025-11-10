# E. coli Genomics Analysis - Complete Example

This example demonstrates a complete E. coli genomics workflow using the WGS Pipeline, based on E. coli Nissle 1917 reference genome analysis.

## Project Setup

### 1. Reference Genome Files

Download E. coli Nissle 1917 reference (GCF_000714595.1):

```bash
# Create reference directory
mkdir -p ~/reference_data/ecoli_nissle1917
cd ~/reference_data/ecoli_nissle1917

# Download from NCBI
datasets download genome accession GCF_000714595.1 --include gff3,protein,genome
unzip ncbi_dataset.zip
cd ncbi_dataset/data/GCF_000714595.1/

# Files needed:
# - GCF_000714595.1_ASM71459v1_genomic.fna (genome)
# - genomic.gff (annotation)
# - protein.faa (proteins)
```

### 2. Data Organization

Organize your sequencing data:

```
~/Data/EcoliProject/
├── WGSData/                          # Main data directory
│   ├── sample1_R1.fastq.gz          # Illumina forward reads
│   ├── sample1_R2.fastq.gz          # Illumina reverse reads
│   ├── sample1_nanopore.fastq.gz    # Nanopore long reads (optional)
│   ├── sample2_R1.fastq.gz
│   ├── sample2_R2.fastq.gz
│   ├── sample2_nanopore.fastq.gz
│   └── samples.txt                   # Sample list file
```

**samples.txt** format (one sample name per line, without file extensions):
```
sample1
sample2
sample3
```

## Complete Analysis Workflow

### Phase 1: Quality Control & Assembly

```bash
# Activate unicycler environment
conda activate unicycler-env

# 1. Nanopore QC (if you have Nanopore reads)
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step nanoplot \
             --threads 32

# 2. Read Trimming
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step trim \
             --threads 32

# 3. FastQC
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step fastqc \
             --threads 32

# 4. Hybrid Assembly (Illumina + Nanopore)
# OR use --step unicycler for Illumina-only
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step hybrid_unicycler \
             --threads 32
```

### Phase 2: Assembly QC & Comparative Genomics

```bash
# Switch to assembly-qc environment
conda activate assembly-qc

# 5. QUAST - Assembly Quality Assessment
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step quast \
             --threads 32

# Expected results for good E. coli assembly:
# - Genome size: ~5.0-5.5 Mb
# - N50: > 1 Mb (ideally complete chromosome)
# - Genome fraction: > 99%
# - Misassemblies: < 5

# 6. Parsnp - Core Genome Alignment & Phylogeny
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step parsnp \
             --threads 32

# Outputs:
# - ~/Data/EcoliProject/Parsnp_Results/parsnp.tree (phylogenetic tree)
# - ~/Data/EcoliProject/Parsnp_Results/parsnp.vcf (core SNPs)
# - ~/Data/EcoliProject/Parsnp_Results/alignment.phylip (alignment)
```

### Phase 3: Genome Annotation

```bash
# Switch to bakta environment
conda activate bakta

# Download Bakta database if not already done
# bakta_db download --output ~/bakta_db --type light

# 7. Bakta Annotation (modern, comprehensive)
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step bakta \
             --bakta_db ~/bakta_db/db \
             --threads 32

# Outputs for each sample:
# - GFF3 files (gene annotations)
# - GenBank files
# - Protein FASTA files
# - TSV tables with functional annotations
```

### Phase 4: Plasmid Analysis

```bash
# Switch to mobsuite environment
conda activate mobsuite

# 8. MOB-suite - Plasmid Identification
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step mobsuite \
             --threads 32

# Results will show:
# - Number of plasmids per strain
# - Plasmid sizes
# - Mobility type (conjugative/mobilizable/non-mobilizable)
# - Incompatibility groups
# - Plasmid FASTA files

# Check summary:
cat ~/Data/EcoliProject/MOBsuite_Results/plasmid_summary.txt
```

### Phase 5: Mutation Analysis (Bacterial Resequencing)

```bash
# Switch to breseq environment
conda activate breseq

# 9. breseq - Mutation Detection
# For Illumina-only analysis:
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step breseq \
             --breseq_mode illumina \
             --threads 32

# For hybrid consensus (Illumina + Nanopore):
# --breseq_mode consensus

# breseq detects:
# - SNPs (point mutations)
# - Small indels
# - Large deletions
# - Mobile element insertions
# - Gene amplifications
# - Mixed populations

# 10. gdtools - Compare Mutations Across Samples
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step gdtools

# View comparison:
# firefox ~/Data/EcoliProject/comparison.html
```

### Phase 6: Traditional Variant Calling (Optional)

```bash
# Switch to wgs-legacy environment
conda activate wgs-legacy

# 11. Snippy - Variant Calling (alternative to breseq)
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step snippy \
             --threads 32

# 12. Snippy-core - Core SNPs
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step snippy_core

# 13. MLST - Sequence Typing
wgs_pipeline --data_dir ~/Data/EcoliProject/WGSData \
             --ref_genome ~/reference_data/ecoli_nissle1917/GCF_000714595.1_ASM71459v1_genomic.fna \
             --ref_gff ~/reference_data/ecoli_nissle1917/genomic.gff \
             --ref_protein ~/reference_data/ecoli_nissle1917/protein.faa \
             --sample_file samples.txt \
             --step mlst
```

## Expected Output Structure

```
~/Data/EcoliProject/
├── WGSData/
│   ├── trimmed_sample1_R1.fastq.gz
│   ├── trimmed_sample1_R2.fastq.gz
│   ├── sample1/
│   │   ├── sample1_NanoPlot/
│   │   ├── sample1_Unicycler_Hybrid/
│   │   └── sample1_FastQC/
│   └── samples.txt
│
├── sample1/
│   ├── sample1_Bakta/           # Bakta annotation
│   ├── sample1_breseq_illumina/ # breseq results
│   └── sample1_Snippy/          # Snippy results
│
├── sample1.fasta                 # Final assembly
├── sample1.gfa                   # Assembly graph
│
├── Parsnp_Results/
│   ├── parsnp.tree              # Phylogenetic tree
│   ├── parsnp.vcf               # Core SNPs
│   └── assemblies/
│
├── MOBsuite_Results/
│   ├── sample1/
│   │   ├── contig_report.txt
│   │   ├── mobtyper_results.txt
│   │   └── plasmid_*.fasta
│   ├── plasmid_fastas/          # All plasmids
│   └── plasmid_summary.txt
│
├── comparison.html               # gdtools comparison
├── GFFFiles/                     # GFF annotations
└── SnippyVCF/                    # Snippy core results
```

## Viewing Results

### 1. Assembly Quality
```bash
# Open QUAST HTML report
firefox ~/Data/EcoliProject/WGSData/../assembly_qc_results/report.html
```

### 2. Phylogenetic Tree
```bash
# View Newick tree
cat ~/Data/EcoliProject/Parsnp_Results/parsnp.tree

# Or use FigTree, iTOL, or ggtree (R)
```

### 3. Mutation Analysis
```bash
# View breseq results for each sample
firefox ~/Data/EcoliProject/sample1/sample1_breseq_illumina/output/index.html

# View multi-sample comparison
firefox ~/Data/EcoliProject/comparison.html
```

### 4. Plasmid Summary
```bash
# Read MOB-suite summary
cat ~/Data/EcoliProject/MOBsuite_Results/plasmid_summary.txt

# View detailed results for a sample
cat ~/Data/EcoliProject/MOBsuite_Results/sample1/contig_report.txt
cat ~/Data/EcoliProject/MOBsuite_Results/sample1/mobtyper_results.txt
```

### 5. Annotations
```bash
# View Bakta annotation table
less ~/Data/EcoliProject/sample1/sample1_Bakta/sample1.tsv

# Search for specific genes (e.g., resistance genes)
grep -i "resistance\|beta-lactam\|ampicillin" ~/Data/EcoliProject/sample1/sample1_Bakta/sample1.tsv
```

## Tips for E. coli Analysis

### Expected Genome Characteristics
- **Chromosome size:** ~4.5-5.5 Mb
- **Plasmids:** 0-5 plasmids, typically 3-200 kb
- **GC content:** ~50-51%
- **Gene count:** ~4,000-5,500 genes
- **Typical N50 for hybrid assembly:** Complete chromosome (5+ Mb)

### When to Use Each Tool

**Parsnp** - Best for:
- Closely related strains (>95% ANI)
- Quick phylogenetic analysis
- Core genome SNP identification

**breseq** - Best for:
- Evolved/passaged strains
- Mixed population detection
- Comprehensive mutation calling
- Gene-level functional impact

**Snippy** - Best for:
- Traditional SNP calling
- Integration with existing workflows
- Quick variant detection

**MOB-suite** - Best for:
- Plasmid epidemiology
- Horizontal gene transfer studies
- Resistance/virulence gene tracking

## Troubleshooting

### Low Assembly Quality
- Check read quality in FastQC
- Verify sufficient coverage (aim for >50x)
- Check for contamination
- Try adjusting Unicycler modes (normal/bold/conservative)

### Parsnp Fails
- Ensure strains are >95% similar
- Check that assemblies are good quality
- Try lowering minimum alignment length
- Verify reference genome is appropriate

### breseq Slow
- Use Illumina-only mode for faster analysis
- Reduce thread count if memory limited
- Consensus mode takes longest but is most accurate

### MOB-suite No Plasmids Found
- Check assembly quality (plasmids may be integrated)
- Verify coverage depth (low coverage may miss plasmids)
- Check contig_report.txt for chromosome vs plasmid classification

## References

- **Unicycler:** Wick RR, et al. (2017) PLoS Comput Biol 13(6): e1005595
- **Parsnp:** Treangen TJ, et al. (2014) Genome Biology 15:524
- **breseq:** Deatherage DE, Barrick JE (2014) Methods Mol Biol 1151:165-188
- **Bakta:** Schwengers et al. (2021) Microb Genom 7(11)
- **MOB-suite:** Robertson & Nash (2018) Microb Genom 4(8)

---

*This example workflow demonstrates the complete analysis pipeline used for E. coli Nissle 1917 variant analysis.*
