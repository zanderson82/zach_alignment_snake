import glob

rule run_cramino:
    input:  
        bam = "".join([PREFIX_REGEX, ".phased.bam"]),
        bai = "".join([PREFIX_REGEX, ".phased.bam.bai"])
    output:
        stats = "".join([PREFIX_REGEX, ".phased.cramino.stats"])
    threads: 10
    conda: "alignment_snake"
    shell:
        """
        cramino -t {threads} {input.bam} > {output.stats}
        """

rule run_samtools_stats:
    input:  
        bam = "".join([WORKDIR, "/", PREFIX_REGEX, ".{type}.bam"]),
        bai = "".join([WORKDIR, "/", PREFIX_REGEX, ".{type}.bam.bai"])
    output:
        stats = temp("".join([WORKDIR, "/", PREFIX_REGEX, ".phased.samtools.stats"]))
    threads: THREADS
    conda: "alignment_snake"
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
        vcf = "".join([PREFIX_REGEX, ".clair3.phased.vcf.gz"]),
        tbi = "".join([PREFIX_REGEX, ".clair3.phased.vcf.gz.tbi"])
    output:
        stats = "".join([PREFIX_REGEX, ".clair3.phased.phasing_stats.tsv"])
    threads: 1
    conda: "alignment_snake"
    shell:
        """
        whatshap stats --tsv={output.stats} {input.vcf}
        """

rule grep_all_cramino:
    input:
        stats = lambda wildcards: glob(("".join(["*-NP-WGS-{project}-*/*-NP-WGS-{project}-*.phased.cramino.stats"])).format(project=config["project"]))
    output:
        out_stats = "".join(["cramino.project.stats.tsv"])
    threads: 1
    shell:
        """
        echo -e "Filename\tYield\tCoverage\tYieldOver25kb\tN50\tIdentity" > {output.out_stats}
        ALLFILES=( {input.stats} )
        for FILE in ${ALLFILES[@]}
        do
            GB=$(sed '4q;d' $FILE | cut -f 2)
            READS=$(sed '2q;d' $FILE | cut -f 2)
            COV=$(sed '5q;d' $FILE | cut -f 2)
            GB25kb=$(sed '6q;d' $FILE | cut -f 2)
            N50=$(sed '7q;d' $FILE | cut -f 2)
            IDENTITY=$(sed '12q;d' $FILE | cut -f 2)
            echo -e "$FILE\t$GB\t$COV\t$GB25kb\t$N50\t$IDENTITY" >> {output.out_stats}
        done
        """

rule grep_phase_stats:
    input:
        stats = lambda wildcards: glob(("".join(["*-NP-WGS-{project}-*/*-NP-WGS-{project}-*.clair3.phased.phasing_stats.tsv"])).format(project=config["project"]))
    output:
        out_stats = "".join(["phasing.project.stats.tsv"])
    threads: 1
    shell:
        """
        ALLFILES=( {input.stats} )
        FILE0 = ${ALLFILES[0]}
        head -n1 $FILE1 > {output.out_stats}
        for FILE in ${ALLFILES[@]}
        do
            grep "ALL" $FILE >> {output.out.stats}
        done
        """