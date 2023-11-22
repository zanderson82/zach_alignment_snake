rule get_target_region_bam:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam.bai"])
    output:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam.bai"])
    threads: int(THREADS/5)
    conda: "alignment_snake"
    params:
        targetBed = lambda wildcards: "".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    shell:
        """
        cat {params.targetBed} | cut -f 2- > tempbed.bed
        samtools view -@ {threads} -M -L tempbed.bed -bo {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        rm tempbed.bed
        """

rule get_region_coverage:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam.bai"])
    output:
        summary = "".join([PREFIX_REGEX, ".phased.target.coverage.tsv"])
    threads: 1
    conda: "rust_plus"
    params:
        targetBed = lambda wildcards: "".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    shell:
        """
        scripts/haplotagStats.sh -i {input.bam} -b {params.targetBed} > {output}
        """

rule run_cramino_target:
    input:  
        bam = "".join([PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([PREFIX_REGEX, ".phased.target.bam.bai"])
    output:
        stats = "".join([PREFIX_REGEX, ".phased.target.cramino.stats"])
    threads: 10
    conda: "alignment_snake"
    shell:
        """
        cramino -t {threads} {input.bam} > {output.stats}
        """