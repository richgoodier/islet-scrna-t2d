#!/bin/bash

# Create output directories
mkdir -p qc_reports trimmed_reads sortmerna_output corrected_reads norm_reads trinity_out spades_out busco_out salmon_quant annotation_results final_transcripts

# Step 1: Quality Control
fastqc *.fastq -o qc_reports/
multiqc qc_reports/ -o qc_reports_summary/

# Step 2: Adapter and Quality Trimming
for file in *.fastq; do
    trim_galore --quality 20 --length 50 --fastqc -o trimmed_reads/ $file
done

# Step 3: Remove rRNA Contamination
for file in trimmed_reads/*.fastq; do
    sortmerna --ref rRNA_database --reads $file --workdir sortmerna_output/
done

# Step 4: Error Correction
for file in sortmerna_output/*.fastq; do
    Rcorrector -t 8 -1 $file -o corrected_reads/
done

# Step 5: Read Normalization
trinityrnaseq/util/insilico_read_normalization.pl --seqType fq --JM 50G --max_cov 50 \
    --left corrected_reads/*_R1.fastq --right corrected_reads/*_R2.fastq --pairs_together --output norm_reads/

# Step 6: Transcriptome Assembly
Trinity --seqType fq --max_memory 100G --left norm_reads/*_R1.fastq --right norm_reads/*_R2.fastq \
    --CPU 16 --output trinity_out/

# Alternative for scRNA-seq
rnaSPAdes.py -1 norm_reads/*_R1.fastq -2 norm_reads/*_R2.fastq -o spades_out/

# Step 7: Assembly Quality Control
TrinityStats.pl trinity_out/Trinity.fasta
busco -m transcriptome -i trinity_out/Trinity.fasta -o busco_out -l eukaryota_odb10 -c 16

# Step 8: Read Mapping and Expression Estimation
salmon index -t trinity_out/Trinity.fasta -i salmon_index
salmon quant -i salmon_index -l A -1 trimmed_reads/*_R1.fastq -2 trimmed_reads/*_R2.fastq \
    -p 8 --validateMappings -o salmon_quant/

# Step 9: Functional Annotation
blastx -query trinity_out/Trinity.fasta -db nr -outfmt 6 -evalue 1e-5 -num_threads 16 -out annotation_results/blast_output.txt
interproscan.sh -i trinity_out/Trinity.fasta -b annotation_results/interpro_out/
eggNOG-mapper -i trinity_out/Trinity.fasta -m diamond -o annotation_results/eggnog_out/

# Step 10: Redundancy Reduction and ORF Prediction
cd-hit-est -i trinity_out/Trinity.fasta -o final_transcripts/non_redundant.fasta -c 0.95
TransDecoder.LongOrfs -t final_transcripts/non_redundant.fasta
TransDecoder.Predict -t final_transcripts/non_redundant.fasta
