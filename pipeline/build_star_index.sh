#!/bin/bash
#SBATCH --job-name=star_index             # Name of the job
#SBATCH --partition=courses               # Who to bill for the job
#SBATCH -N 1                              # How many nodes do you need. When in doubt, use 1
#SBATCH -c 16                             # How many "threads" do you need for your job
#SBATCH --mem 32G                         # How much memory
#SBATCH -t 8:00:00                        # How long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email>
#SBATCH --out=./genome_index/logs/%x_%j.log
#SBATCH --error=./genome_index/logs/%x_%j.err

# Setup environment
conda activate islet-scrna-t2d

# Create logs directory if it doesn't exist
mkdir -p ./genome_index/logs

# Run STAR index generation
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir ./genome_index \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.109.gtf \
     --sjdbOverhang 100