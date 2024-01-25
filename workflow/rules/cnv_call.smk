rule run_qdnaseq:
    input:
        aligned_unphased_bam="".join([SAMPLE_WORKPATH, ".notPhased.bam"]),
        aligned_bam_index = "".join([SAMPLE_WORKPATH, ".notPhased.bam.bai"])
    output:
        qdnaseq_seg=temp("".join([SAMPLE_WORKPATH, ".called_cnv.seg"]))
        qdnaseq_bins=temp("".join([SAMPLE_WORKPATH, ".cnv.bins.txt"]))
        qdnaseq_vcf="".join([SAMPLE_WORKPATH, ".called_cnv.vcf"])
        qdnaseq_pdf="".join([SAMPLE_WORKPATH, ".called_cnv.pdf"])
    threads: THREADS
    conda:
         config["conda_r"]
    log:
        o = "".join(["logs/",LOG_REGEX,"run_qdnaseq","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"run_qdnaseq","-stderr.log"])
    params:
        output_trunk=SAMPLE_WORKPATH
    shell:
        """
       RScript workflow/scripts/qdnaseq_sop.R {input.aligned_unphased_bam} {params.output_trunk} {threads} > {log.o} 2> {log.e}
        """

rule plot_qdnaseq:
    input:
        qdnaseq_seg="".join([SAMPLE_WORKPATH, ".called_cnv.seg"])
        qdnaseq_bins="".join([SAMPLE_WORKPATH, ".cnv.bins.txt"])
    output:
        qdnaseq_plot="".join([SAMPLE_WORKPATH, ".called_cnv.detail_plot.pdf"])
    conda:
        config["conda_r"]
    log:
        o = "".join(["logs/",LOG_REGEX,"plot_qdnaseq","-stdout.log"]),
        e = "".join(["logs/",LOG_REGEX,"plot_qdnaseq","-stderr.log"])
    shell:
        """
        RScript workflow/scripts/qdnaseq_plot.R {input.qdnaseq_seg} {input.qdnaseq_bins} {output.qdnaseq_plot}
        """