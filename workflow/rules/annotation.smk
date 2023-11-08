# this rule should run vep and then make the output more readable.
rule run_vep:
    input:
        clair3_phased_vcf="".join([SAMPLE_WORKPATH, ".clair3.phased.vcf.gz"])
    output:
        vep_vcf=temp("".join([SAMPLE_WORKPATH, ".clair3.phased.vep.vcf"])),
    threads: THREADS
    log:
        o = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-stdout.log",
        e = "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-stderr.log"
    params:
        CADD="/data/dat/annotationData/CADDv1.6_hg38_whole_genome_SNVs.tsv.gz",
        SPLICEAISNV="/data/dat/annotationData/spliceai_scores.raw.snv.hg38.vcf.gz",
        SPLICEAIINDEL="/data/dat/annotationData/spliceai_scores.raw.indel.hg38.vcf.gz",
        GNOMAD="/data/dat/annotationData/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
        CLINVAR="/data/dat/annotationData/clinvar.hg38.vcf.gz",
        cache_directory="/home/dm1/.vep",
        plugin_dir="/home/dm1/.vep/Plugins"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        vep --cache --merged --offline --dir_cache {params.cache_directory} --force_overwrite --vcf --pick --use_given_ref --fork {THREADS} -i {input.clair3_phased_vcf} -o {output.vep_vcf} --canonical --symbol --numbers --domains --pubmed --sift b --polyphen b --regulatory --total_length --numbers --af --max_af --af_1kg --dir_plugins {params.plugin_dir} --plugin CADD,{params.CADD} --plugin SpliceAI,snv={params.SPLICEAISNV},indel={params.SPLICEAIINDEL} --custom {params.GNOMAD},gnomADg,vcf,exact,0,AF --custom {params.CLINVAR},ClinVar,vcf,exact,0,CLINSIG,CLNREVSTAT,CLNDN 2>> {log.e}
        """

# this rule  should filter the vep vcf
rule filter_vep:
    input:
        vep_vcf="".join([SAMPLE_WORKPATH, ".clair3.phased.vep.vcf"])
    output:
        vep_intermediate=temp("".join([SAMPLE_WORKPATH, ".vep.vcf.tmp"])),
        vep_lt1_phased_vcf=temp("".join([SAMPLE_WORKPATH, ".clair3.phased.vep.af_lt_1_phased.csv"])),
        vep_lt1_notPhased_vcf=temp("".join([SAMPLE_WORKPATH, ".clair3.phased.vep.af_lt_1_notPhased.csv"]))
        
    log:
        "logs/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}-{MB}-stdout.log"
    conda:
         "../envs/alignment.yaml"
    script:
        "../scripts/filter_vep.py"
