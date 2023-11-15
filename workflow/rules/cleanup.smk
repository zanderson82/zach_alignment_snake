rule move_vcf:
    input:
        vcf = "".join([WORKDIR, "/", PREFIX_REGEX, "{caller}.{phasing}.vcf"])
    output:
        vcf = protected("".join([PREFIX_REGEX, "{caller}.{phasing}.vcf"]))
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        """

rule move_vcf_gz:
    input:
        vcf = "".join([WORKDIR, "/", PREFIX_REGEX, "{caller}.{phasing}.vcf.gz"]),
        tbi = "".join([WORKDIR, "/", PREFIX_REGEX, "{caller}.{phasing}.vcf.gz.tbi"])
    output:
        vcf = protected("".join([PREFIX_REGEX, "{caller}.{phasing}.vcf.gz"])),
        tbi = protected("".join([PREFIX_REGEX, "{caller}.{phasing}.vcf.gz.tbi"]))
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.tbi} {output.tbi}
        """

rule move_bam:
    input:
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".{type}.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".{type}.bam.bai"])
    output:
        bam = protected("".join([PREFIX_REGEX, ".{type}.bam"])),
        bai = protected("".join([PREFIX_REGEX, ".{type}.bam.bai"]))
    threads: 1
    shell:
        """
        cp {input.bam} {output.bam}
        cp {input.bai} {output.bai}
        """


rule move_csv:
    input:
        vep_phased = "".join([WORKDIR, "/", PREFIX_REGEX, ".clair3.phased.vep.af_lt_1_phased.csv"]),
        vep_notphased = "".join([WORKDIR, "/",  PREFIX_REGEX, ".clair3.phased.vep.af_lt_1_notPhased.csv"])
    output:
        vep_phased = protected("".join([PREFIX_REGEX, ".clair3.phased.vep.af_lt_1_phased.csv"])),
        vep_notphased = protected("".join([PREFIX_REGEX, ".clair3.phased.vep.af_lt_1_notPhased.csv"]))
    threads: 1
    shell:
        """
        cp {input.vep_phased} {output.vep_phased}
        cp {input.vep_notphased} {output.vep_notphased}
        """