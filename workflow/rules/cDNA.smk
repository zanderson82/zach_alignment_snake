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
    threads: THREADS
    conda:
         config["conda_alignment"]
    log:
        o = "".join(["logs/",LOG_REGEX,"pychop","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"pychop","-stderr.log"])
    params:
        kit=get_kit
    shell:
        """
        echo "minimap2 -ax splice -t {THREADS} -y {REFGENOME} {input.full_read_set} | samtools sort -@ {THREADS}" >> {log.o}
        minimap2 -ax splice -t {THREADS} -y {REFGENOME} {input.full_read_set} | samtools sort -@ {THREADS} > {output.aligned_bam} 2> {log.e}
        """

# stringtie

rule runStringTie:
    input:
        aligned_bam = "".join([SAMPLE_WORKPATH, ".pychop.aligned.bam"])
    output:
        aligned_bam = temp("".join([SAMPLE_WORKPATH, ".pychop.aligned.bam"])),
    threads: THREADS
    conda:
         config["conda_stringtie"]
    log:
        o = "".join(["logs/",LOG_REGEX,"pychop","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"pychop","-stderr.log"])
    params:
        kit=get_kit
    shell:
        """
        echo "minimap2 -ax splice -t {THREADS} -y {REFGENOME} {input.full_read_set} | samtools sort -@ {THREADS}" >> {log.o}
        minimap2 -ax splice -t {THREADS} -y {REFGENOME} {input.full_read_set} | samtools sort -@ {THREADS} > {output.aligned_bam} 2> {log.e}
        """

# gffcompare

# gffread



