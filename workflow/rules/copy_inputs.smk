rule copy_alignedBam_to_workdir:
    input: 
        bam="".join([PREFIX_REGEX,".phased.bam"]),
        index="".join([PREFIX_REGEX,".phased.bam.bai"])
    output: 
        bam=temp("".join([WORKDIR,"/",PREFIX_REGEX,".phased.bam"])),
        index=temp("".join([WORKDIR,"/",PREFIX_REGEX,".phased.bam.bai"]))
    threads: 1
    shell:
        """
        cp {input.bam} {output.bam}
        cp {input.index} {output.index}
        """

rule copy_clairVCF_to_workdir:
    input: 
        vcf="".join([PREFIX_REGEX,".clair3.{phasing}.vcf.gz"]),
        index="".join([PREFIX_REGEX,".clair3.{phasing}.vcf.gz.tbi"])
    output: 
        vcf=temp("".join([WORKDIR,"/",PREFIX_REGEX,".clair3.{phasing}.vcf.gz"])),
        index=temp("".join([WORKDIR,"/",PREFIX_REGEX,".clair3.{phasing}.vcf.gz.tbi"]))
    threads: 1
    shell:
        """
        cp {input.vcf} {output.vcf}
        cp {input.index} {output.index}
        """
