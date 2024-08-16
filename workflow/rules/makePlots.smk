## Read Lengths plot

rule subsample_bam:
    input: "".join([PREFIX_REGEX, ".phased.bam"])
    output: 
        bam=temp("".join([PREFIX_REGEX, ".subsampled.phased.bam"])),
        bai=temp("".join([PREFIX_REGEX, ".subsampled.phased.bam.bai"])),
        stats=temp("".join([PREFIX_REGEX, ".subsampled.phased.stats"]))
    threads: THREADS
    params:
        f=0.10
    conda: config["conda_samtools"]
    shell:
        """
        samtools view -@ {THREADS} --subsample {params.f} -bo {output.bam} {input}  
        samtools index -@ {THREADS} {output.bam}
        samtools stats -@ {THREADS} {output.bam} > {output.stats}
        """

rule plot_readLengths:
    input: 
        stats="".join([PREFIX_REGEX, ".subsampled.phased.stats"])
    output: 
        df=temp("".join([PREFIX_REGEX, ".subsampled.phased.readlengths.tsv"])),
        plot=temp("".join([PREFIX_REGEX, ".plot_readlengths.png"]))
    threads: 1
    conda: config["conda_r"]
    params:
        color="cornflowerblue",
        script="workflow/scripts/plotSamtoolsLength.R"
    shell:
        """
        echo -e "Length\tCount" > {output.df}
        grep ^RL {input.stats} | cut -f2,3 >> {output.df}

        Rscript {params.script} {output.df} {output.plot} {params.color}
        """

rule plot_ru_readLengths:
    input: 
        stats="".join([PREFIX_REGEX, ".phased.target.samtools.stats"])
    output: 
        df=temp("".join([PREFIX_REGEX, ".phased.target.readlengths.tsv"])),
        plot=temp("".join([PREFIX_REGEX, ".target.plot_readlengths.png"]))
    threads: 1
    conda: config["conda_r"]
    params:
        color="darkorchid3",
        script="workflow/scripts/plotSamtoolsLength.R"
    shell:
        """
        echo -e "Length\tCount" > {output.df}
        grep ^RL {input.stats} | cut -f2,3 >> {output.df}

        Rscript {params.script} {output.df} {output.plot} {params.color}
        """

## Coverage Plot

rule plot_coverage:
    input: "".join([PREFIX_REGEX, ".windowed.coverage.txt"])
    output: temp("".join([PREFIX_REGEX, ".plot_depth_coverage.png"]))
    threads: 1
    params: 
        positions="",
        script="workflow/scripts/plotCoverage.R"
    conda: config["conda_karyoplotR"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule get_windowed_coverage:
    input: 
        "".join([PREFIX_REGEX, ".phased.bam"])
    output: 
        temp("".join([PREFIX_REGEX, ".windowed.coverage.txt"]))
    threads: 1
    params: 
        positions="/n/dat/hg38/hg38.500kb.windowed.positions"
    conda: config["conda_samtools"]
    shell:
        """
        samtools mpileup -q 1 --positions {params.positions} --output-extra MAPQ {input} | cut -f1,2,4,7 > {output}
        """

## SNP and INDEL quality plots

rule prepare_ClairtoPlot:
    input: "".join([PREFIX_REGEX, ".clair3.notPhased.vcf.gz"])
    output: temp("".join([PREFIX_REGEX, ".clair3.notPhased.forPlotting.vcf"]))
    threads: 1
    shell:
        """
        INPUT={input}
        tempinput=${{INPUT%*.vcf.gz}}
        cp $INPUT $tempinput.forPlotting.vcf.gz

        gunzip $tempinput.forPlotting.vcf.gz
        """

rule make_clair_qual_plots:
    input: "".join([PREFIX_REGEX, ".clair3.notPhased.forPlotting.vcf"])
    output: 
        plotdf=temp("".join([PREFIX_REGEX, ".clair3.notPhased.forQualityPlotting.tsv"])),
        indelplot=temp("".join([PREFIX_REGEX, ".plot_indel_quality.png"])),
        snpplot=temp("".join([PREFIX_REGEX, ".plot_snv_quality.png"]))
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


# SNP and INDEL quality plots for targeted region

rule prepare_ru_ClairtoPlot:
    input: "".join([PREFIX_REGEX, ".clair3.notPhased.vcf.gz"])
    output: 
        vcf=temp("".join([PREFIX_REGEX, ".target.clair3.notPhased.forPlotting.vcf"])),
        tempbed=temp("".join([PREFIX_REGEX, ".nameless.bed"]))
    threads: 1
    conda: config["conda_bcftools"]
    params:
        targetBed = lambda wildcards: "".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    shell:
        """
        cut -f 2,3,4 {params.targetBed} > {output.tempbed}
        bcftools view -R {output.tempbed} {input} > {output.vcf}
        """

rule make_clair_qual_target_plots:
    input: "".join([PREFIX_REGEX, ".target.clair3.notPhased.forPlotting.vcf"])
    output: 
        plotdf=temp("".join([PREFIX_REGEX, ".target.clair3.notPhased.forQualityPlotting.tsv"])),
        indelplot=temp("".join([PREFIX_REGEX, ".target.plot_indel_quality.png"])),
        snpplot=temp("".join([PREFIX_REGEX, ".target.plot_snv_quality.png"]))
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

rule make_efficiency_plots:
    input: "".join([PREFIX_REGEX, ".phased.cramino.stats"])
    output: 
        plotdf=temp("".join([PREFIX_REGEX, ".efficiencyStats.tsv"])),
        readplot=temp("".join([PREFIX_REGEX, ".read_efficiency.png"])),
        yieldplot=temp("".join([PREFIX_REGEX, ".yield_efficiency.png"])),
        n50plot=temp("".join([PREFIX_REGEX, ".n50_efficiency.png"]))
    threads: 1
    conda: config["conda_r"]
    params:
        df_script="workflow/scripts/efficiencyPlot.sh",
        plot_script="workflow/scripts/plotEfficiency.R",
        basecall_path="/data/prealign_qc/dorado_summary/qual_only",
        rawdf=temp("".join([PREFIX_REGEX, ".rawSeq.tsv"])),
        libraryfolder=get_basecall_folder
    shell:
        """
        bash {params.df_script} -c {input} -l '{params.libraryfolder}' -o {output.plotdf} -r {params.rawdf} -b {params.basecall_path}
        Rscript {params.plot_script} {output.plotdf} {output.readplot} Reads
        Rscript {params.plot_script} {output.plotdf} {output.yieldplot} Yield
        Rscript {params.plot_script} {output.plotdf} {output.n50plot} N50
        """


#### Phasing plots

## SNP plot

rule setup_clair_snp_plot:
    input: "".join([PREFIX_REGEX, ".clair3.notPhased.forPlotting.vcf"])
    output: 
        plotdf=temp("".join([PREFIX_REGEX, ".clair3.notPhased.forSNPlotting.tsv"])),
        bindf=temp("".join([PREFIX_REGEX, ".clair3.notPhased.forSNPlotting.bins.tsv"])),
        filterdf=temp("".join([PREFIX_REGEX, ".clair3.notPhased.forSNPlotting.1MBfilter.tsv"])),
        vcf=temp("".join([PREFIX_REGEX, ".clair3.notPhased.lowHeterozygosity.vcf"]))
    threads: 1
    params:
        plot_script="workflow/scripts/plotSNPs.R",
        ru_plot_script="workflow/scripts/plotSNPsTargets.R",
        THRESHOLD=0.35,
        COUNTTHRESH=50,
        MINGTQUAL=10
    shell:
        """
        echo -e "Chromosome\tPosition\tGenotype\tGTQuality" > {output.plotdf}
        grep -v ^# {input} | cut -f1,2,10 | tr ':' '\t' | cut -f1,2,3,4 | awk -v mingt={params.MINGTQUAL} '{{if($4 > mingt){{print $0}}}}' | tr ' ' '\t' >> {output.plotdf}

        echo "creating bins"

        awk 'NR==2{{chr=$1;pos=$2;if($3=="1/1"){{\
        hom=1;het=0}}\
        else{{\
        hom=0;het=1}};\
        binstart=$2}}\
        NR>2{{\
        if($2+0 < (binstart+10000) && $1==chr){{\
            if($3=="1/1"){{\
            hom+=1}}\
            else{{\
            het+=1\
            }};\
            pos=$2\
        }}else{{\
        print chr,binstart,pos,hom/(hom+het);\
        chr=$1;pos=$2;binstart=$2;\
        if($3=="1/1"){{\
            hom=1\
            }}else{{\
            het=1}}\
        }}\
        }}\
        END{{print chr,binstart,pos,hom/(hom+het)}}' {output.plotdf} | tr ' ' ',' > {output.bindf}

        echo "iterating over bins"

        echo -e "Chromosome\tStart\tStop\tThreshold\tBins\tLength\tAverageHomozygosity" > {output.filterdf}
        cat {output.bindf} | tr ',' '\t' | awk -v thresh={params.THRESHOLD} 'NR==1{{chr=$1;\
        start=$2;\
        stop=$3;\
        if($4 > thresh){{toggle="above"}}else{{toggle="below"}};\
        count=1;\
        total=$4}}\
        NR>1{{if($1 != chr){{print chr,start,stop,toggle,count,stop-start,total/count;\
        chr=$1;start=$2;stop=$3;count=1;total=$4;\
        if($4>thresh){{toggle="above"}}else{{toggle="below"}}}}\
        else{{\
        if($4 >= thresh && toggle=="above"){{\
        count+=1;stop=$3;total+=$4}}\
        if($4 >= thresh && toggle=="below"){{\
        print chr,start,stop,toggle,count,stop-start,total/count;\
        toggle="above";start=$2;stop=$3;count=1;total=$4;}}\
        if($4 < thresh && toggle=="below"){{\
        count+=1;stop=$3;total+=$4;}}\
        if($4 < thresh && toggle=="above"){{\
        print chr,start,stop,toggle,count,stop-start,total/count;\
        toggle="below";start=$2;stop=$3;count=1;total=$4}}\
        }}}}\
        END{{print chr,start,stop,toggle,count,stop-start,total/count}}' | awk -v thresh={params.COUNTTHRESH} '{{if($4=="above" && $5 >= thresh){{print $0}}}}' | \
        awk 'NR==1{{chr=$1;start=$2;stop=$3;bins=$5;hom=$7;}}\
        NR>1{{if(chr!=$1){{print chr,start,stop,"above",bins,stop-start,hom;\
        chr=$1;start=$2;stop=$3;bins=$5;hom=$7}}else{{\
        if($2-5000000 < stop){{stop=$3;hom=(($7*$5)+(bins*hom))/(bins+$5);bins+=$5}}else{{\
        print chr,start,stop,"above",bins,stop-start,hom;\
        chr=$1;start=$2;stop=$3;bins=$5;hom=$7}}\
        }}\
        }}\
        END{{print chr,start,stop,"above",bins,stop-start,hom}}' | tr ' ' '\t' >> {output.filterdf}

        echo "setting up VCF"
        echo "##fileformat=VCFv4.2" > {output.vcf}
        echo "##source=/n/scripts/clairAnalysis.sh" >> {output.vcf}
        grep "^##reference" {input} >> {output.vcf}
        grep "^##contig" {input} | grep -v chrUn | grep -v random | grep -v chrM >> {output.vcf}
        echo '##REF=<ID=DIP,Description="Low heterozygosity call">' >> {output.vcf}
        echo '##ALT=<ID=HOM,Description="Homozygous stretch, proportion > 0.4 length > 100 bins">'  >> {output.vcf}
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of variant">' >> {output.vcf}
        echo '##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of variant">' >> {output.vcf}
        echo '##INFO=<ID=FREQ,Number=1,Type=Float,Description="Average proportion of homozygous SNVs">' >> {output.vcf}
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{input}" >> {output.vcf}

        echo "populating VCF"

        cut -f1,2,3,5,6,7 {output.filterdf} | awk '{{print $1,$2,".","<DIP>","<HOM>",$4,"PASS","END="$3";LEN="$5";FREQ="$6,"GT","1/1"}}' | tr ' ' '\t' >> {output.vcf}

        """

rule make_clair_snp_plot:
    input:
        plotdf="".join([PREFIX_REGEX, ".clair3.notPhased.forSNPlotting.tsv"]),
        filterdf="".join([PREFIX_REGEX, ".clair3.notPhased.forSNPlotting.1MBfilter.tsv"])
    output: 
        plot=temp("".join([PREFIX_REGEX, ".clair3.notPhased.snpPlot.png"]))
    threads: 1
    params:
        multiallelic=0,
        targetBed = get_target_bed
    conda: config["conda_r"]
    shell:
        """
        strategy={wildcards.STRATEGY}
        if [ $strategy == "RU" ]
        then
            regionstring=$( awk '{{print $2":"$3"-"$4}}' {params.targetBed} | tr '\n' ';')
            Rscript workflow/scripts/plotSNPsTargets.R {input.plotdf} {output} {params.multiallelic} {input.filterdf} $regionstring
        else
            Rscript workflow/scripts/plotSNPs.R {input.plotdf} {output} {params.multiallelic} {input.filterdf}
        fi
        """

## karyotype with het/hom phasing

rule plot_whatshap:
    input: "".join([PREFIX_REGEX, ".clair3.phased.phasing_stats.tsv"])
    output: plot=temp("".join([PREFIX_REGEX, ".clair3.phased.whatshap_plot.png"]))
    threads: 1
    conda: config["conda_r"]
    params: 
        script=" workflow/scripts/plotWhatshap.R"
    shell:
        """
        Rscript {params.script} {input} {output.plot}
        """


## haplotagging and phaseblock membership

rule run_hp_dp_long:
    input:  
        bam = "".join([PREFIX_REGEX, ".phased.bam"]),
        bai = "".join([PREFIX_REGEX, ".phased.bam.bai"])
    output:
        stats = temp("".join([PREFIX_REGEX, ".longform.hp_dp.stats"])),
        temp_positions=temp("".join([PREFIX_REGEX, ".temp_hpdp.positions"])),
        temp_pileup=temp("".join([PREFIX_REGEX, ".temp_hpdp.pileup.tsv"]))
    threads: 1
    conda: config["conda_rust"]
    params: 
        targets = get_target_bed,
        temp_prefix="".join([PREFIX_REGEX, ".temp_hpdp"])
    shell:
        """
        bash workflow/scripts/haplotagStats.sh -i {input.bam} -b {params.targets} -l -t {params.temp_prefix} > {output.stats}
        """


rule plot_hpdp_detail:
    input: "".join([PREFIX_REGEX, ".longform.hp_dp.stats"])
    output: "".join([PREFIX_REGEX, ".hp_dp_long_complete.txt"])
    threads: 1
    conda: config["conda_karyoplotR"]
    params: 
        script=" workflow/scripts/plotHPDP.R",
        output_trunk="".join([PREFIX_REGEX, ".hp_dp.detail_plot.png"])
    shell:
        """
        Rscript {params.script} {input} {params.output_trunk} {output}
        """
rule generate_hpdp_table:
    input: "".join([PREFIX_REGEX, ".hp_dp.stats"])
    output: 
        html = temp("".join([PREFIX_REGEX, ".hp_dp.snippet.html"]))
    threads: 1
    shell:
        """
        bash workflow/scripts/csvToTable.sh  -i {input} -o {output.html} -f "2,3,4,5,6,7,8"
        """

## SV counts

## VEP summary

rule filter_vep_table:
    input: "".join([PREFIX_REGEX, ".clair3.phased.vep.111.af_lt_1.csv"])
    output: temp("".join([PREFIX_REGEX, ".clair3.phased.vep.111.af_lt_1.pathogenic.csv"]))
    threads: 1
    shell:
        """
        grep -i pathogenic {input} | grep -vi "benign" | awk -F ',' '{{if($77>20){{print$0}}}}' | awk -F ',' '{{if($41 != "."){{print $0}}}}' > {output}
        """

rule filter_vep_target_table:
    input: "".join([PREFIX_REGEX, ".clair3.phased.vep.111.af_lt_1.csv"])
    output: temp("".join([PREFIX_REGEX, ".clair3.phased.vep.111.af_lt_1.target.csv"]))
    params: targetGenes=get_target_genes
    threads: 1
    shell:
        """
        head -n1 {input} > {output}
        GENES=( {params.targetGenes} )
        for GENE in ${{GENES[@]}}
        do
            grep ",$GENE", {input} | grep -vi "benign" | sort -t ',' -k77n,77 >> {output}
        done
        """


rule generate_vep_table:
    input: "".join([PREFIX_REGEX, ".clair3.phased.vep.111.af_lt_1.{type}.csv"])
    output: 
        html = temp("".join([PREFIX_REGEX, ".vep.{type}.snippet.html"]))
    threads: 1
    shell:
        """
        bash workflow/scripts/csvToTable.sh  -i {input} -o {output.html} -f "1,2,7,8,9,10,11,12,13,15,17,26,28,29,42,43,51,73,75,77,89,91,93"
        """

## Fill report

rule make_report:
    input: get_report_inputs
    output: "".join([PREFIX_REGEX, ".alignment_report.html"])
    threads: 1
    params:
        html_template = branch( evaluate("{STRATEGY}=='RU'"), then="workflow/resources/ru_template_report.html", otherwise="workflow/resources/report_template.html"),
        script = "workflow/scripts/fill_report.sh",
        bamfiles = get_target_bams,
        libname=PREFIX,
        server=config["server"],
        email=config["email"]
    shell:
        """
        cp {params.html_template} {output}
        BAMARRAY=( {params.bamfiles} )
        BAMINPUT=$( echo ${{BAMARRAY[@]}} | tr ' ' ',' )
        CRAMINO={params.libname}/{params.libname}.phased.cramino.stats
        HPDP={params.libname}/{params.libname}.hp_dp.stats
        if [ {wildcards.STRATEGY} == "RU" ]
        then
            bash {params.script} -o {output} -l {params.libname} -b $BAMINPUT -c $CRAMINO -h $HPDP -r
        else
            bash {params.script} -o {output} -l {params.libname} -b $BAMINPUT -c $CRAMINO -h $HPDP
        fi
        file={output}
        filename=${{file##*/}}
        if [ {params.server} == "franklin" ]
        then
            echo "Library {params.libname} has been aligned and a report is attached." | mail --content-name $filename -A {output} -s "Alignment complete: {params.libname}" {params.email}
        else
            echo "Library {params.libname} has been aligned and a report is attached." | mail -a {output} -s "Alignment complete: {params.libname}" {params.email}
        fi
        """