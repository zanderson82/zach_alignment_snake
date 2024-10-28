rule get_target_region_bam:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam.bai"])
    output:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam.bai"]),
        tempbed = temp("".join([WORKDIR, "/", PREFIX_REGEX, ".tempbed"]))
    threads: int(THREADS/5)
    conda: config["conda_samtools"]
    params:
        targetBed = lambda wildcards: "".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    shell:
        """
        cat {params.targetBed} | cut -f 2- > {output.tempbed}
        samtools view -@ {threads} -M -L {output.tempbed} -bo {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule get_region_coverage:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.target.bam.bai"])
    output:
        summary = "".join([PREFIX_REGEX, ".phased.target.coverage.tsv"])
    threads: 1
    conda: config["conda_rust"]
    params:
        targetBed = lambda wildcards: "".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    shell:
        """
        workflow/scripts/haplotagStats.sh -i {input.bam} -b {params.targetBed} > {output}
        """

rule run_cramino_target:
    input:  
        bam= f"{FINALDIR}/{PREFIX_REGEX}.phased.target.bam",
        bai= f"{FINALDIR}/{PREFIX_REGEX}.phased.target.bam.bai"
    output:
        stats = f"{FINALDIR}/{PREFIX_REGEX}.phased.target.cramino.stats"
    threads: 10
    conda: config["conda_cramino"]
    shell:
        """
        cramino -t {threads} {input.bam} > {output.stats}
        """

rule run_samtools_target:
    input:  
        bam= f"{FINALDIR}/{PREFIX_REGEX}.phased.target.bam",
        bai= f"{FINALDIR}/{PREFIX_REGEX}.phased.target.bam.bai"
    output:
        stats = f"{FINALDIR}/{PREFIX_REGEX}.phased.target.samtools.stats"
    threads: 10
    conda: config["conda_samtools"]
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output.stats}
        """

rule run_hp_dp_target:
    input:  
        bam= f"{FINALDIR}/{PREFIX_REGEX}.phased.target.bam",
        bai= f"{FINALDIR}/{PREFIX_REGEX}.phased.target.bam.bai"
    output:
        stats = f"{FINALDIR}/{PREFIX_REGEX}.target.hp_dp.stats"
    threads: 1
    conda: config["conda_rust"]
    params:
                targets=get_hpdp_png_names
    shell:
        """
        bash workflow/scripts/haplotagStats.sh -i {input.bam} -b {params.targets} > {output.stats}
        """