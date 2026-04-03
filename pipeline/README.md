# Pipeline — HPC Scripts

This directory contains the HPC pipeline for processing raw FASTQ files through alignment and read counting. These steps require a cluster environment (tested on SLURM).

## Prerequisites

Install into your conda environment (see root `environment.yml`), plus:

```bash
conda install -c bioconda star cutadapt nextflow subread
```

Reference genome files (download separately):

- `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- `Homo_sapiens.GRCh38.109.gtf`

## Steps

### 1. Build STAR genome index

Requires ~32 GB RAM and ~2 hours. Submit as a SLURM job:

```bash
sbatch build_star_index.sh
```

Edit `build_star_index.sh` to set `GENOME_DIR`, `FASTA`, and `GTF` paths for your cluster.

### 2. Trim adapters

```bash
nextflow run fastq_processing.nf
```

Runs cutadapt (quality threshold 20) in parallel on all FASTQ files.

### 3. Align to genome

```bash
nextflow run star_align.nf
```

Aligns trimmed reads to GRCh38 using STAR (project used v2.3.0e). Processed 149 FASTQ files in ~3.5 hours on our cluster.

### 4. Download scripts

`shell_scripts/` contains the ENA download scripts used to retrieve raw FASTQ files. See [`../samples/README.md`](../samples/README.md) for sample selection details.

## Output

Aligned BAM files are passed to featureCounts (see [`../gene_matrix/README.md`](../gene_matrix/README.md)).
