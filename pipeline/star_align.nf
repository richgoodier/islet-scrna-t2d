nextflow.enable.dsl = 2

params.reads = '/scratch/goodier.r/project_files/trimmed_t2d/*.fastq'
// params.reads = '/scratch/goodier.r/project_files/trimmed_healthy/*.fastq'

params.outdir = '/scratch/goodier.r/project_files/trimmed_t2d/aligned_t2d'
// params.outdir = '/scratch/goodier.r/project_files/trimmed_healthy/aligned_healthy'

params.genomeDir = '/scratch/goodier.r/project_files/genome/genome_index'

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
            --outTmpDir ${sample_id}_STARtmp \
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