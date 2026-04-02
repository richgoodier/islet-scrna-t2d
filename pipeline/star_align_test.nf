nextflow.enable.dsl = 2

params.reads = '/courses/BINF6310.202530/students/SEC04_Group_3/files/test/trimmed_outputs/*.fastq'
params.genomeDir = '/courses/BINF6310.202530/students/SEC04_Group_3/files/test/genome/genome_index1'
params.outdir = '/courses/BINF6310.202530/students/SEC04_Group_3/files/test/aligned'
params.threads = 8

Channel.fromPath(params.reads)
    .set { trimmed_fastqs }

process STAR_Align {
    input:
        tuple path(fastq), val(genomeDir)

    output:
        path "*.bam"

    script:
        sample_id = fastq.getBaseName().replaceFirst(/^trimmed_/, "").replaceFirst(/\.fastq$/, "")
        """
        STAR --genomeDir $genomeDir \
            --readFilesIn $fastq \
            --runThreadN ${params.threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${sample_id}_

        # Rename STAR output to something simpler
        mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam

        # Copy results to final output dir
        mkdir -p ${params.outdir}
        cp ${sample_id}.bam ${params.outdir}/
        
        cp -f ${sample_id}.Log.final.out ${params.outdir}/ || true
        """
}

workflow {
    trimmed_fastqs
        .combine(Channel.value(params.genomeDir))
        | STAR_Align
}