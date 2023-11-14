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
        o = "".join(["logs/",LOG_REGEX,"-{svcaller}-","phase_vcf","-stdout.log"])
        e = "".join(["logs/",LOG_REGEX,"-{svcaller}-","phase_vcf","-stderr.log"])
    params:
        OUTPUT_DIR=get_output_dir
    conda:
         "../envs/alignment.yaml"
    shell:
        """
        longphase_tmp={params.OUTPUT_DIR}/longphase_tmp_{wildcards.svcaller}
        longphase phase -s {input.clair3_unphased_vcf} -b {input.aligned_unphased_bam} -r {REFGENOME} -t {THREADS} --sv-file={input.sv_vcf} -o $longphase_tmp --ont 2>> {log.e}
        mv "$longphase_tmp"_SV.vcf {output.sv_phased_vcf}
        rm "$longphase_tmp".vcf
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
         "../envs/alignment.yaml"
    threads: THREADS
    log:
        o = "".join(["logs/",LOG_REGEX,"phase_bamfile","-stdout.log"])
        e = "".join(["logs/",LOG_REGEX,"phase_bamfile","-stderr.log"])
    params:
        qualityThreshold=1
    shell:
        """
        longphase haplotag --snp-file={input.clair3_phased_vcf} --bam-file={input.aligned_unphased_bam} --qualityThreshold={params.qualityThreshold} -t {THREADS} --sv-file={input.sniffles_phased} -o {output.phased_bam} 2>> {log.e}
        mv {output.phased_bam}.bam {output.phased_bam}
        samtools index -@ {THREADS} {output.phased_bam} 2>> {log.e}
        """