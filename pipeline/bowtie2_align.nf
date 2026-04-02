nextflow.enable.dsl = 2
// Define the input channel
inputFile = Channel.fromPath('./trimmed_*.fastq')

// Define the alignment process
process align {
input:
path fastq
output:
path "alignment.sam"
script:
"""
bowtie2 -x ${params.reference_index} -U ${fastq} -S alignment.sam
"""
}

// Define pipeline parameters
params.reference_index = '/courses/BINF6310.202530/'

// Define the workflow
workflow {
inputFile | align
}
