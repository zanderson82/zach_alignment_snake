# pychopper
def get_kit(wildcards):
    return samples.loc[wildcards.SAMPLEID, "Kit"]

rule pychop:
    input:
        fastqfile = "".join([SAMPLE_WORKPATH, ".fastq"]),
        fastq_completion="".join([SAMPLE_WORKPATH, "-temp_fastq.log"])
    output:
        full_read_set = temp("".join([SAMPLE_WORKPATH, ".pychop.fastq"]))
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

rule skipPychop:
    input:
        fastqfile = "".join([SAMPLE_WORKPATH, ".fastq"]),
        fastq_completion="".join([SAMPLE_WORKPATH, "-temp_fastq.log"])
    output:
        full_read_set = temp("".join([SAMPLE_WORKPATH, ".raw.fastq"]))
    threads: THREADS
    shell:
        """
        mv {input.fastqfile} {output.full_read_set}
        """

# align

rule alignCDNA:
    input:
        full_read_set = "".join([SAMPLE_WORKPATH, ".{chopped}.fastq"])
    output:
        aligned_bam = temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.bam"])),
        aligned_index = temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.bam.bai"]))
    threads: THREADS
    conda:
         config["conda_alignment"]
    log:
        o = "".join(["logs/",LOG_REGEX,".{chopped}.aligncDNA","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,".{chopped}.aligncDNA","-stderr.log"])
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
        aligned_bam = "".join([SAMPLE_WORKPATH, ".{chopped}.aligned.bam"])
    output:
        gtf = temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gtf"])),
        abundance = temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.abundance.tab"])),
        coverage = temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.coverage.gtf"]))
    threads: THREADS
    conda:
         config["conda_stringtie"]
    log:
        o = "".join(["logs/",LOG_REGEX,".{chopped}.stringtie","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,".{chopped}.stringtie","-stderr.log"])
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
        gtf="".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gtf"])
    output:
        annotation=temp( "".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.annotated.gtf"]) ),
        refmap=temp( "".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.refmap"]) ),
        bestrefmap=temp( "".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.tmap"]) ),
        tracking=temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.tracking"])),
        loci=temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.loci"])),
        log=temp("".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare"]))
    threads: THREADS
    conda:
        config["conda_gffcompare"]
    log:
        o = "".join(["logs/",LOG_REGEX,".{chopped}.gffcompare","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,".{chopped}.gffcompare","-stderr.log"])
    params:
        transcript_ref=config["transcript_reference"],
        prefix="".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare"])
    shell:
        """
        gffcompare -V -R -r {params.transcript_ref} -o {params.prefix} {input.gtf}
        infile="{input.gtf}"
        filename="${{infile##*/}}"
        mv {params.prefix}.$filename.refmap {output.refmap}
        mv {params.prefix}.$filename.tmap {output.bestrefmap}
        #cat gffcompare >> {log.o}
        """


# gffread

rule gffRead:
    input:
        annotated_gtf="".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.annotated.gtf"])
    output:
        transcriptome=temp( "".join([SAMPLE_WORKPATH, ".{chopped}.aligned.gffcompare.annotated.transcriptome.fa"]) )
    threads: THREADS
    conda:
        config["conda_gffread"]
    log:
        o = "".join(["logs/",LOG_REGEX,".{chopped}.gffread","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,".{chopped}.gffread","-stderr.log"])
    params:
        genome_ref=config["refgenome"]
    shell:
        """
        gffread -w {output.transcriptome} -g {params.genome_ref} 
        """

