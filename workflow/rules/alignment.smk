rule make_fastqs:
    input:
        get_target_bams
    output:
        fastqfile = temp("".join([SAMPLE_WORKPATH, ".fastq"])),
        fastq_completion=temp("".join([SAMPLE_WORKPATH, "-temp_fastq.log"]))
    log:
        o = "".join(["logs/",LOG_REGEX,"make_fastqs","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"make_fastqs","-stderr.log"])
    threads: THREADS
    conda:
        config["conda_samtools"]
    params:
        tags = "MM,ML",
        outdir=get_output_dir,
        finaldir=get_final_dir,
        email=config["email"]
    shell:
        """
        echo "Alignment for {wildcards.SAMPLEID} in project 1000genomes is starting. FASTQ files are being generated from {input}" | mail -s  "Alignment for {wildcards.SAMPLEID} Starting" {params.email}
        mkdir -p {params.outdir}
        mkdir -p {params.finaldir}
        inbams=( {input} )
        for bam in ${{inbams[@]}}
        do
            echo "running samtools fastq with $bam" >> {log.e}
            samtools fastq --threads {threads} -T {params.tags} $bam >> {output.fastqfile} 2>> {log.e}
        done
        touch {output.fastq_completion}
        """

rule make_alignment:
    input:
        fastqfile = "".join([SAMPLE_WORKPATH, ".fastq"]),
        fastq_completion="".join([SAMPLE_WORKPATH, "-temp_fastq.log"])
    output:
        aligned_bam = temp("".join([SAMPLE_WORKPATH, ".notPhased.bam"])),
    threads: THREADS
    conda:
         config["conda_minimap"]
    log:
        o = "".join(["logs/",LOG_REGEX,"make_alignment","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"make_alignment","-stderr.log"])
    shell:
        """
        echo running 'minimap2 -t {THREADS} -y -a {ONTMMIFILE} {input.fastqfile} 2>> {log.e} | samtools sort -@ {THREADS} -o {output.aligned_bam} 2>> {log.e}' >> {log.o}
        minimap2 -t {THREADS} -y -a {ONTMMIFILE} {input.fastqfile} 2>> {log.e} | samtools sort -@ {THREADS} -o {output.aligned_bam} 2>> {log.e}
        """

rule index_alignment:
    input:
        aligned_bam = "".join([SAMPLE_WORKPATH, ".notPhased.bam"])
    output:
        aligned_bam_index = temp("".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"]))
    threads: THREADS
    conda:
         config["conda_samtools"]
    shell:
        """
        samtools index -@ {THREADS} {input.aligned_bam}
        """

rule run_clair3:
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"])
    output:
        clair3_phased_vcf=temp("".join([SAMPLE_WORKPATH, ".clair3.phased.vcf.gz"])),
        clair3_phased_vcf_index=temp("".join([SAMPLE_WORKPATH, ".clair3.phased.vcf.gz.tbi"])),
        clair3_not_phased_vcf=temp("".join([SAMPLE_WORKPATH, ".clair3.notPhased.vcf.gz"])),
        clair3_not_phased_vcf_index=temp("".join([SAMPLE_WORKPATH, ".clair3.notPhased.vcf.gz.tbi"]))
    threads: THREADS
    log:
        o = "".join(["logs/",LOG_REGEX,"run_clair3","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"run_clair3","-stderr.log"])
    params:
        OUTPUT_DIR=get_output_dir,
        cmodel=get_clair_model,
        chromosomes=config["chromosomes"]
    conda:
         config["conda_clair3"]
    shell:
        """echo "running clair3" >> {log.o}
        run_clair3.sh --ctg_name={params.chromosomes} --bam_fn={input.aligned_unphased_bam} --ref_fn={REFGENOME} --threads={THREADS} --platform=ont --model_path={params.cmodel} --output={params.OUTPUT_DIR} --enable_phasing 2>> {log.e}
        mv {params.OUTPUT_DIR}/merge_output.vcf.gz {output.clair3_not_phased_vcf}
        mv {params.OUTPUT_DIR}/merge_output.vcf.gz.tbi {output.clair3_not_phased_vcf_index}
        mv {params.OUTPUT_DIR}/phased_merge_output.vcf.gz {output.clair3_phased_vcf}
        mv {params.OUTPUT_DIR}/phased_merge_output.vcf.gz.tbi {output.clair3_phased_vcf_index}
        rm -rf {params.OUTPUT_DIR}/tmp
        """