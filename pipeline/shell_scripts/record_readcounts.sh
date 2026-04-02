#!/bin/bash

# Output header
echo "directory,filename,read_count" > fastq_readcount.csv

# Find all .fastq files recursively
find . -type f -name "*.fastq" | while read file; do
    read_count=$(wc -l < "$file")
    read_count=$((read_count / 4))
    directory=$(dirname "$file")
    filename=$(basename "$file")
    echo "$directory,$filename,$read_count" >> fastq_readcount.csv
done
