rule get_target_region_bam:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam.bai"])
    output:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam.bai"])
    threads: THREADS
    conda: "../envs/alignment.yaml"
    params:
        targetBed = lambda wildcards: samples.loc[wildcards.SAMPLEID, "BedFile"]
    shell:
        """
        samtools view -@ {threads} -M -L {params.targetBed} {threads} {input.bam} > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule get_region_coverage:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam.bai"])
    output:
        summary = "".join([PREFIX_REGEX, ".phased.target.coverage.tsv"])
    threads: 1
    conda: "../envs/rustenv.yaml"
    params:
        targetBed = lambda wildcards: samples.loc[wildcards.SAMPLEID, "BedFile"]
    shell:
        """
        ../scripts/haplotagStats.sh -i {input.bam} -b {params.targetBed} > {output}
        """

rule run_cramino_target:
    input:  
        bam = "".join([PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([PREFIX_REGEX, ".phased.target.bam.bai"])
    output:
        stats = "".join([PREFIX_REGEX, ".phased.target.cramino.stats"])
    threads: 10
    conda: "../envs/alignment.yaml"
    shell:
        """
        cramino -t {threads} {input.bam} > {output.stats}
        """