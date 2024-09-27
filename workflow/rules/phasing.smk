# this rule uses Longphase to phase structural variant vcfs.
rule phase_vcf:
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"]),
        clair3_unphased_vcf="".join([SAMPLE_WORKPATH, ".clair3.notPhased.vcf.gz"]),
        sv_vcf="".join([SAMPLE_WORKPATH, ".sv_{svcaller}.notPhased.vcf"])
    output:
        sv_phased_vcf=temp("".join([SAMPLE_WORKPATH, ".sv_{svcaller}.phased.vcf"]))
    threads: THREADS
    log:
        o = "".join(["logs/",LOG_REGEX,"-{svcaller}-","phase_vcf","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"-{svcaller}-","phase_vcf","-stderr.log"])
    params:
        OUTPUT_DIR=get_output_dir
    conda:
        config["conda_longphase"]
    shell:
        """
        longphase_tmp={params.OUTPUT_DIR}/longphase_tmp_{wildcards.svcaller}
        longphase phase -s {input.clair3_unphased_vcf} -b {input.aligned_unphased_bam} -r {REFGENOME} -t {THREADS} --sv-file={input.sv_vcf} -o $longphase_tmp --ont 2>> {log.e}
        mv "$longphase_tmp"_SV.vcf {output.sv_phased_vcf}
        rm "$longphase_tmp".vcf
        """

rule phase_indel_vcf:
    input:
        vcf="".join([SAMPLE_WORKPATH, ".clair3.notPhased.vcf.gz"]),
        bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        bai="".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"])
    output:
        vcf=temp("".join([SAMPLE_WORKPATH, ".whatshap.phased_indels.vcf.gz"]))
    threads: THREADS
    log:
        o = "".join(["logs/",LOG_REGEX,"-whatshap_rephase-","phase_vcf","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"-whatshap_rephase-","phase_vcf","-stderr.log"])
    conda:
        config["conda_clair3"]
    shell:
        """
        whatshap phase -o {output.vcf} --reference={REFGENOME} --ignore-read-groups --indels {input.vcf} {input.bam}
        """

rule index_indel_vcf:
    input:
        vcf="".join([SAMPLE_WORKPATH, ".whatshap.phased_indels.vcf.gz"]),
    output:
        index=temp("".join([SAMPLE_WORKPATH, ".whatshap.phased_indels.vcf.gz.csi"]))
    threads: THREADS
    log:
        o = "".join(["logs/",LOG_REGEX,"-whatshap_rephase-","phase_vcf","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"-whatshap_rephase-","phase_vcf","-stderr.log"])
    conda:
        config["conda_bcftools"]
    shell:
        """
        bcftools index {input.vcf}
        """

# this rule uses the phased sniffles sv vcf, the phased clair snp vcf, and Longphase to phase the unphased bam file, then indexes it.
rule phase_bamfile:
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"]),
        sniffles_phased="".join([SAMPLE_WORKPATH, ".sv_sniffles.phased.vcf"]),
        clair3_phased_vcf="".join([SAMPLE_WORKPATH, ".clair3.phased.vcf.gz"])
    output:
        phased_bam=temp("".join([SAMPLE_WORKPATH, ".phased.bam"])),
        phased_bam_index=temp("".join([SAMPLE_WORKPATH, ".phased.bam.bai"]))
    conda:
        config["conda_longphase"]
    threads: THREADS
    log:
        o = "".join(["logs/",LOG_REGEX,"phase_bamfile","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"phase_bamfile","-stderr.log"])
    params:
        qualityThreshold=1,
        refgenome=REFGENOME
    shell:
        """
        longphase haplotag --snp-file={input.clair3_phased_vcf} --bam-file={input.aligned_unphased_bam} --qualityThreshold={params.qualityThreshold} -t {THREADS} --sv-file={input.sniffles_phased} -o {output.phased_bam} -r {params.refgenome} 2>> {log.e}
        mv {output.phased_bam}.bam {output.phased_bam}
        samtools index -@ {THREADS} {output.phased_bam} 2>> {log.e}
        """