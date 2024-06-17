## Read Lengths plot
rule subsample_bam:
    input: "".join([SAMPLE_WORKPATH, ".phased.bam"])
    output: 
        bam=temp("".join([SAMPLE_WORKPATH, ".subsampled.phased.bam"])),
        bai=temp("".join([SAMPLE_WORKPATH, ".subsampled.phased.bam.bai"])),
        stats=temp("".join([SAMPLE_WORKPATH, ".subsampled.phased.stats"]))
    threads: THREADS
    params:
        f=0.10
    conda: config["conda_alignment"]
    shell:
        """
        samtools view -@ {THREADS} --subsample {params.f} -bo {output.bam} {input}  
        samtools index {output.bam}
        samtools stats -@ {THREADS} {output.bam} > {output.stats}
        """

rule plot_readLengths:
    input: 
        stats="".join([SAMPLE_WORKPATH, ".subsampled.phased.stats"])
    output: 
        df=temp("".join([SAMPLE_WORKPATH, ".subsampled.phased.readlengths.tsv"])),
        plot=temp("".join([SAMPLE_WORKPATH, ".plot_readlengths.png"]))
    threads: 1
    conda: config["conda_r"]
    params:
        color="cornflowerblue",
        script="workflow/scripts/plotSamtoolsLength.R"
    shell:
        """
        echo -e "Length\tCount" > {output.df}
        grep ^RL {input.stats} | cut -f2,3 >> {output.df}

        Rscript {params.script} {input.stats} {output.plot} {params.color}
        """

## Coverage Plot

rule get_windowed_coverage:
    input: "".join([SAMPLE_WORKPATH, ".phased.bam"])
    output: temp("".join([SAMPLE_WORKPATH, ".windowed.coverage.txt"]))
    threads: 1
    params: 
        positions=""
    conda: config["conda_alignment"]
    shell:
        """
        samtools mpileup -q 1 --positions {params.positions} --output-extra MAPQ {input} | cut -f1,2,4,7 > {output}
        """

rule plot_coverage:
    input: "".join([SAMPLE_WORKPATH, ".windowed.coverage.txt"])
    output: temp("".join([SAMPLE_WORKPATH, ".plot_depth_coverage.txt"]))
    threads: 1
    params: 
        positions=""
    conda: config["conda_R"]
    shell:
        """
        """

## SNP and INDEL quality plots

rule prepare_ClairtoPlot:
    input: "".join([SAMPLE_WORKPATH, ".clair3.notPhased.vcf.gz"])
    output: temp("".join([SAMPLE_WORKPATH, ".clair3.notPhased.forPlotting.vcf"]))
    threads: 1
    shell:
        """
        INPUT={input}
        tempinput=${{INPUT%*.vcf.gz}}
        cp $INPUT $tempinput.forPlotting.vcf.gz

        gunzip $tempinput.forPlotting.vcf.gz
        """

rule make_clair_qual_plots:
    input: "".join([SAMPLE_WORKPATH, ".clair3.notPhased.forPlotting.vcf"])
    output: 
        plotdf=temp("".join([SAMPLE_WORKPATH, ".clair3.notPhased.forQualityPlotting.tsv"]))
        indelplot=temp("".join([SAMPLE_WORKPATH, ".plot_indel_quality.png"]))
        snpplot=temp("".join([SAMPLE_WORKPATH, ".plot_snv_quality.png"]))
    threads: 1
    conda: config["conda_r"]
    params:
        plot_script="workflow/scripts/plotClairQuality.R"
    shell:
        """
        echo -e "Chromosome\tPosition\tType\tQuality" > {output.plotdf}
        grep -v ^# {input} | cut -f1,2,4,5,6 | awk '{{if($4 ~ /,/ ){{next}}else{{if( length($4) != length($5)){{print $1,$2,"INDEL",$5}}else{{\
        print $1,$2,"SNV",$5}}}}}}' | tr ' ' '\t' >> {output.plotdf}

        Rscript {params.plot_script} {output.plotdf} {output.indelplot} INDEL darkseagreen
        Rscript {params.plot_script} {output.plotdf} {output.snpplot} SNV darksalmon
        """


###### Efficiency Plots

## Reads


## GB

## N50

#### Phasing plots

## SNP plot

## karyotype with het/hom phasing

## haplotagging and phaseblock membership
