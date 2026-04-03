import os
import gzip
import csv

def count_reads_in_fastq(file_path):
    """Counts the number of reads in a FASTQ or FASTQ.GZ file."""
    open_func = gzip.open if file_path.endswith(".gz") else open
    with open_func(file_path, 'rt') as f:
        line_count = sum(1 for _ in f)
    return line_count // 4  # Each read consists of 4 lines

def process_directory(directory, output_csv):
    """Processes all .fastq and .fastq.gz files in a directory and saves read counts to a CSV file."""
    results = []
    
    for file in os.listdir(directory):
        if file.endswith(".fastq") or file.endswith(".fastq.gz"):
            file_path = os.path.join(directory, file)
            read_count = count_reads_in_fastq(file_path)
            results.append([file, read_count])
            print(f"Processed {file}: {read_count} reads")

    # Save results to a CSV file
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Filename", "Read Count"])
        writer.writerows(results)

    print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    directory = r"/path/to/fastq/files"
    output_csv = "fastq_read_counts.csv"
    process_directory(directory, output_csv)
