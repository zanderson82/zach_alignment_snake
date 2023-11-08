# add clean up or clean up at end
rule run_sniffles2:
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"])
    output:
        sniffles_unphased=temp("".join([SAMPLE_WORKPATH, ".sv_sniffles.notPhased.vcf"]))
    threads: THREADS
    conda:
         "../envs/alignment.yaml"
    log:
        o = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-sniffles-stdout.log",
        e = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-sniffles-stderr.log"
    shell:
        """
        sniffles --input {input.aligned_unphased_bam} --output-rnames --vcf {output.sniffles_unphased} --threads {THREADS} --allow-overwrite 2>> {log.e}
        """


rule run_svim: 
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"])
    output:
        svim_unphased=temp("".join([SAMPLE_WORKPATH, ".sv_svim.notPhased.vcf"]))
    threads: 1
    log:
        o = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-svim-stdout.log",
        e = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-svim-stderr.log"
    params:
        OUTPUT_DIR=get_output_dir
    conda:
         "../envs/alignment.yaml"
    shell:
        """
        outSvim="{params.OUTPUT_DIR}/svim_output"
        [ -d $outSvim ] && rm -r $outSvim
        svim alignment --read_names $outSvim {input.aligned_unphased_bam} {REFGENOME} 2>> {log.e}
        cat $outSvim/variants.vcf | grep \"^#\" | grep -v chrUn | grep -v random > {output.svim_unphased} 2>> {log.e}
        cat $outSvim/variants.vcf | grep -v \"^#\" | grep -v chrUn | grep -v random >> {output.svim_unphased} 2>> {log.e}
        """


# add clean up or clean up at end
rule run_cuteSV:
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"])
    output:
        cutesv_unphased=temp("".join([SAMPLE_WORKPATH, ".sv_cutesv.notPhased.vcf"]))
    threads: THREADS
    log:
        o = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-cutesv-stdout.log",
        e = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-cutesv-stderr.log"
    params:
        OUTPUT_DIR=get_output_dir
    conda:
         "../envs/alignment.yaml"
    shell:
        """
        outCuteSV="{params.OUTPUT_DIR}/cuteSV_work_dir"
        [ -d $outCuteSV ] && rm -r $outCuteSV
        mkdir $outCuteSV
        cuteSV {input.aligned_unphased_bam} {REFGENOME} {output.cutesv_unphased} "$outCuteSV" --threads {THREADS} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --report_readid --genotype 2>> {log.e}
        """
        
