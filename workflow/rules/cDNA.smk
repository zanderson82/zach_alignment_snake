# pychopper
def get_kit(wildcards):
    return samples.loc[SAMPLEID, "kit"]

rule pychop:
    input:
        fastqfile = "".join([SAMPLE_WORKPATH, ".fastq"]),
        fastq_completion="".join([SAMPLE_WORKPATH, "-temp_fastq.log"])
    output:
        full_read_set = temp("".join([SAMPLE_WORKPATH, ".pychop.bam"]))
    threads: THREADS
    conda:
         config["conda_pychopper"]
    log:
        o = "".join(["logs/",LOG_REGEX,"pychop","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"pychop","-stderr.log"])
    params:
        kit=get_kit
    shell:
        """
        echo "running pychopper -k {params.kit} -r -t {THREADS} -U -y {input.fastqfile} {output.full_read_set}" >> {log.o}
        pychopper -k {params.kit} -r -t {THREADS} -U -y {input.fastqfile} {output.full_read_set} 2> {log.e}
        """

# align

rule alignCDNA:
    input:
        full_read_set = "".join([SAMPLE_WORKPATH, ".pychop.bam"])
    output:
        aligned_bam = temp("".join([SAMPLE_WORKPATH, ".pychop.aligned.bam"]))
        aligned_bam = temp("".join([SAMPLE_WORKPATH, ".pychop.aligned.bam.bai"]))
    threads: THREADS
    conda:
         config["conda_alignment"]
    log:
        o = "".join(["logs/",LOG_REGEX,"aligncDNA","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"aligncDNA","-stderr.log"])
    params:
        kit=get_kit
    shell:
        """
        echo "minimap2 -ax splice -t {THREADS} -y {REFGENOME} {input.full_read_set} | samtools sort -@ {THREADS}" >> {log.o}
        minimap2 -ax splice -t {THREADS} -y {REFGENOME} {input.full_read_set} | samtools sort -@ {THREADS} > {output.aligned_bam} 2> {log.e}
        samtools index -@ {THREADS} {output.aligned_bam}
        """

# stringtie

rule runStringTie:
    input:
        aligned_bam = "".join([SAMPLE_WORKPATH, ".pychop.aligned.bam"])
    output:
        gtf = temp("".join([SAMPLE_WORKPATH, ".pychop.aligned.gtf"])),
        abundance = temp("".join([SAMPLE_WORKPATH, ".pychop.aligned.abundance.tab"])),
        coverage = temp("".join([SAMPLE_WORKPATH, ".pychop.aligned.coverage.gtf"])),
    threads: THREADS
    conda:
         config["conda_stringtie"]
    log:
        o = "".join(["logs/",LOG_REGEX,"stringtie","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"stringtie","-stderr.log"])
    params:
        transcript_ref=config["transcript_reference"],
        genome_ref=config["refgenome"]
    shell:
        """
        echo "stringtie -L -o {output.gtf} -G {params.transcript_ref} -p {THREADS} -A {output.abundance} -C {output.coverage} --ref {params.genome_ref} {input.aligned_bam}" >> {log.o}
        stringtie -L -o {output.gtf} -G {params.transcript_ref} -p {THREADS} -A {output.abundance} -C {output.coverage} --ref {params.genome_ref} {input.aligned_bam} 2> {log.e}
        """

# gffcompare

rule gffCompare:
    input:
        gtf="".join([SAMPLE_WORKPATH, ".pychop.aligned.gtf"])
    output:
        stats=temp( "".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare.stats"]) ),
        annotation=temp( "".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare.annotated.gtf"]) ),
        refmap=temp( "".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare.refmap"]) ),
        bestrefmap=temp( "".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare.tmap"]) )
    threads: THREADS
    conda:
        config["conda_gffcompare"]
    log:
        o = "".join(["logs/",LOG_REGEX,"gffcompare","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"gffcompare","-stderr.log"])
    params:
        transcript_ref=config["transcript_reference"]
        prefix="".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare"])
    shell:
        """
        gffcompare -V -R -r {params.transcript_ref} -o {params.prefix} {input.gtf}
        infile="{input.gtf}"
        filename="${{infile##*/}}"
        mv {params.prefix}.$filename.refmap {output.refmap}
        mv {params.prefix}.$filename.tmap {output.bestrefmap}
        """


# gffread

rule gffRead:
    input:
        annotated_gtf="".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare.annotated.gtf"])
    output:
        transcriptome=temp( "".join([SAMPLE_WORKPATH, ".pychop.aligned.gffcompare.annotated.transcriptome.fa"]) )
    threads: THREADS
    conda:
        config["conda_gffread"]
    log:
        o = "".join(["logs/",LOG_REGEX,"gffread","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"gffread","-stderr.log"])
    params:
        genome_ref=config["refgenome"]
    shell:
        """
        gffread -w {output.transcriptome} -g {params.genome_ref} 
        """

