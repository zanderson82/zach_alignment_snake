import glob

rule run_cramino:
    input:  
        bam = f"{FINALDIR}/{PREFIX_REGEX}.phased.bam",
        bai = f"{FINALDIR}/{PREFIX_REGEX}.phased.bam.bai"
    output:
        stats = f"{FINALDIR}/{PREFIX_REGEX}.phased.cramino.stats"
    threads: 10
    conda: config["conda_cramino"]
    shell:
        """
        cramino -t {threads} {input.bam} > {output.stats}
        """

rule run_hp_dp:
    input:  
        bam = f"{FINALDIR}/{PREFIX_REGEX}.phased.bam",
        bai = f"{FINALDIR}/{PREFIX_REGEX}.phased.bam.bai"
    output:
        stats = f"{FINALDIR}/{PREFIX_REGEX}.hp_dp.stats"
    threads: 1
    conda: config["conda_rust"]
    params: 
        targets=get_target_bed
    shell:
        """
        bash workflow/scripts/haplotagStats.sh -i {input.bam} -b {params.targets} > {output.stats}
        """
        

rule run_samtools_stats:
    input:  
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.bam.bai"])
    output:
        stats = temp("".join([WORKDIR, "/", PREFIX_REGEX, ".phased.samtools.stats"]))
    threads: THREADS
    conda: config["conda_samtools"]
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output.stats}
        """

rule grep_samtools_stats:
    input:  
        stats = "".join([WORKDIR, "/", PREFIX_REGEX, ".phased.samtools.stats"])
    output:
        summary = "".join([PREFIX_REGEX, ".phased.samtools.stat.summary"])
    threads: 1
    shell:
        """
        INPUTFILE={input.stats}
        TOTBASES=$(grep "total length:" $INPUTFILE | cut -f 3)
        TOTGB=$(echo "scale=2; $TOTBASES/1000000000" | bc)
        TOTREADS=$(grep "SN\tsequences:" $INPUTFILE | cut)
        MQ0=$(grep "SN\treads MQ0:" $INPUTFILE | cut)
        MAPPINGERROR=$(grep "SN\terror rate:" $INPUTFILE | cut)
        MEANQUAL=$(grep "SN\taverage quality:" $INPUTFILE | cut)
        MEAN_RL=$(grep "SN\taverage length:" $INPUTFILE | cut)
        ##
        MAPQs=$( grep "^MAPQ" $INPUTFILE | awk 'BEGIN{total=0}{total+=($2*$3)}END{print total}')
        MIDQ=$(echo "$MAPQs /2" | bc)
        MAPQ_MEDIAN=$( grep "^MAPQ" $INPUTFILE | awk -v tl="$MIDQ" 'BEGIN{{sum=0}}{{sum+=$3;if(sum>tl){{print $2; exit 0;}}}}')

        MID=$(echo "$TOTBASES /2" | bc)
        N50=$( awk '/^RL/ {{for(i=0;i<$3;i++){{print $2}}}}' "$INPUTFILE" | awk -v tl="$MID" 'NR==1{{sum=$1}}NR>1{{sum=sum+$1;if(sum>tl){{print $1; exit 0;}}}}')
        N50kb=$(echo "scale=2; $N50/1000" | bc)
        """

rule run_whatshap:
    input:
        vcf = f"{FINALDIR}/{PREFIX_REGEX}.clair3.phased.vcf.gz",
        tbi = f"{FINALDIR}/{PREFIX_REGEX}.clair3.phased.vcf.gz.tbi"
    output:
        stats = f"{FINALDIR}/{PREFIX_REGEX}.clair3.phased.phasing_stats.tsv"
    threads: 1
    conda: config["conda_clair3"]
    shell:
        """
        whatshap stats --tsv={output.stats} {input.vcf}
        """
