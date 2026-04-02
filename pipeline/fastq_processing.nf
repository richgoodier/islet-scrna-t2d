nextflow.enable.dsl = 2

params.quality = 20

// Define the input channel
fastqFiles = Channel.fromPath('./*.fastq')


// Define the trimming process
process trimAndFilter {
    input:
    path fastq

    output:
    path "trimmed/trimmed_${fastq.baseName}.fastq"

    script:
    """
    mkdir -p trimmed
    cutadapt -q ${params.quality} \\
             -o trimmed/trimmed_${fastq.baseName}.fastq \\
             ${fastq} > trimmed/${fastq.baseName}_cutadapt.log
    """
}

// Define the workflow
workflow {
    fastqFiles | trimAndFilter
}