import pandas as pd
import subprocess as sp
import glob

THREADS=config["threads"]
REFGENOME=config["refgenome"]
ONTMMIFILE=config["ontmmifile"]
WORKDIR=config["working_dir"]
FINALDIR=config["final_dir"]
INDIR=config["input_dir"]
PROJECT=config["project"]
BEDFILEDIR=config["bedfiledir"]

basecalled_bam_string=config["basecalled_bam_string"]

PREFIX=config["prefix_regex"]
lookupPrefix=config["prefix_lookup"]
prefix_column_format=config["sampleDB_prefix_format_columns"].split(",")
prefix_column_names=config["sampleDB_prefix_column_names"].split(",")
index_col=prefix_column_names[0]

samples = pd.read_table(config["samples"], sep="\t").set_index(index_col)

PREFIX_REGEX="/".join([PREFIX,PREFIX])
SAMPLE_WORKPATH="".join([WORKDIR, "/", PREFIX_REGEX]) # SAMPLE_WORKPATH should be used in all inputs up to moving the files to /n/alignments
LOG_REGEX=PREFIX

def apply_suffix(wildcards, ending, sampleID):
    # applies a suffix to the path of a sample folder and sample prefix, e.g. M1079-NP-{STRATEGY}-Project-OutsideID-pb/M1079-NP-{STRATEGY}-Project-OutsideID-pb.vcf.gz
    return ".".join([get_output_name(wildcards, sampleID), ending])

def get_prefix(wildcards, sampleID):
    if len(prefix_column_format) != len(prefix_column_names):
        raise ValueError("sampleDB_prefix_format_columns must contain the same number of comma separated values as sampleDB_prefix_column_names")
    sampleVals=[sampleID]
    sampleVals+=samples.loc[sampleID][prefix_column_names[1:]].values.tolist()
    column_pairs=dict(zip(prefix_column_format, sampleVals))
    return lookupPrefix.format(**column_pairs)

def get_output_name(wildcards, sampleID):
    #returns the path to a sample folder and sample prefix, e.g. M1079-NP-{STRATEGY}-Project-OutsideID-pb/M1079-NP-{STRATEGY}-ProjectOutsideID-pb
    prefix=get_prefix(wildcards, sampleID)
    return "/".join([prefix, prefix])

def get_output_dir(wildcards):
    # returns the working (Franklin) directory with specific folder, an absolute path. e.g. /data/alignments/M1079-NP-{STRATEGY}-Project-OutsideID-pb/
    sampleID=wildcards.SAMPLEID
    prefix=get_prefix(wildcards, sampleID)
    return "/".join(["{outdir}", prefix]).format(outdir=WORKDIR)

def get_final_dir(wildcards):
    # returns the destination (McClintock) directory with specific folder, a relative path
    sampleID=wildcards.SAMPLEID
    prefix=get_prefix(wildcards, sampleID)
    return "/".join(["{finaldir}", prefix]).format(finaldir=FINALDIR)

def get_report_targets(wildcards):
    f = open(config["targetfile"], "r")
    targets = f.read().split("\n")
    f.close()
    final_targets=[]
    if config["explicitLibraries"]:
        targetsamples=[x.split("-")[0] for x in targets]
    else:
        targetsamples=targets
    final_targets = [apply_suffix(wildcards, "alignment_report.html", ts) for ts in targetsamples]
    return final_targets

def get_targets_new(wildcards):
    f = open(config["targetfile"], "r")
    targets = f.read().split("\n")
    f.close()
    final_targets=[]
    if config["explicitLibraries"]:
        targetsamples=[x.split("-")[0] for x in targets]
    else:
        targetsamples=targets
    final_targets = []
    if config["qcCaller"] != "cramino" and config["qcCaller"] != "samtools":
        print("qcCaller value {} not recognized, using cramino".format(config["qcCaller"]))
        summarizer="cramino"
    else:
        summarizer=config["qcCaller"]
    endings=[]
    if config["outputs"]["alignBam"] or config["allTargets"]:
        endings+=["phased.bam", "phased.bam.bai"]
    if config["outputs"]["clair3"] or config["allTargets"]:
        endings+=["clair3.phased.vcf.gz", "clair3.notPhased.vcf.gz", "clair3.phased.vcf.gz.tbi", "clair3.notPhased.vcf.gz.tbi"]
    if config["outputs"]["sniffles"] or config["allTargets"]:
        endings+=["sv_sniffles.phased.vcf", "sv_sniffles.notPhased.vcf"]
    if config["outputs"]["svim"] or config["allTargets"]:
        endings+=["sv_svim.phased.vcf", "sv_svim.notPhased.vcf"]
    if config["outputs"]["cuteSV"] or config["allTargets"]:
        endings+=["sv_cutesv.phased.vcf", "sv_cutesv.notPhased.vcf"]
    if config["outputs"]["CNVcalls"] or config["allTargets"]:
        endings+=["called_cnv.vcf", "called_cnv.pdf", "called_cnv.detail_plot.png","called_cnv.annotated.vcf"]
    if config["outputs"]["VEP"] or config["allTargets"]:
        endings+=["clair3.phased.vep.111.vcf", "clair3.phased.vep.111.af_lt_1.csv"]
    if config["outputs"]["basicQC"] or config["allTargets"]:
        endings+=["phased.{}.stats".format(summarizer)]
    if config["outputs"]["phaseQC"] or config["allTargets"]:
        endings+=["clair3.phased.phasing_stats.tsv"]
    if config["outputs"]["report"] or config["allTargets"]:
        endings+=["alignment_report.html"]
    for ts in targetsamples:
        strategy=samples.loc[ts,"Strategy"]
        file_endings=endings
        if strategy == "RU":
            file_endings+=["phased.target.bam", "phased.target.bam.bai"]
            if config["outputs"]["basicQC"] or config["allTargets"]:
                file_endings+=["phased.target.{}.stats".format(summarizer)]
            if config["outputs"]["phaseQC"] or config["allTargets"]:
                file_endings+=["hp_dp.stats"]
        else:
            file_endings+=["hp_dp.stats"]
        all_targets = [apply_suffix(wildcards, x, ts) for x in file_endings]
        final_targets += all_targets
    return final_targets
        

def get_targets_transcriptome(wildcards):
    f = open(config["targetfile"], "r")
    targets = f.read().split("\n")
    f.close()
    final_targets=[]
    if config["explicitLibraries"]:
        targetsamples=[x.split("-")[0] for x in targets]
    else:
        targetsamples=targets
    final_targets = []
    endings=[]
    if config["transcriptomeOutputs"]["alignment"] or config["allTranscriptomeTargets"]:
        endings+=["aligned.bam", "aligned.bam.bai"]
    if config["transcriptomeOutputs"]["annotation"] or config["allTranscriptomeTargets"]:
        endings+=["aligned.abundance.tab", "aligned.coverage.gtf", "aligned.gffcompare.annotated.gtf", "aligned.gffcompare.refmap", "aligned.gffcompare.tmap"]
    if config["transcriptomeOutputs"]["transcriptome"] or config["allTranscriptomeTargets"]:
        endings+=["aligned.gffcompare.annotated.transcriptome.fa"]
    if config["transcriptomeOutputs"]["pychop"] or config["allTranscriptomeTargets"]:
        endings=["pychop."+x for x in endings]
    else:
        endings=["raw."+x for x in endings]
    for ts in targetsamples:
        file_endings=endings
        all_targets = [apply_suffix(wildcards, x, ts) for x in file_endings]
        final_targets += all_targets
    return final_targets
        
def get_target_bams(wildcards):
    if config["explicitLibraries"]:
        f = open(config["targetfile"], "r")
        targets = f.read().split("\n")
        f.close()
        libraries=list(filter(lambda x: x.split("-")[0]==wildcards.SAMPLEID, targets))
        return ["{}/{}".format(INDIR,x) for x in libraries]
    else:
        folder="/".join([INDIR,basecalled_bam_string.format(wildcards=wildcards)])
        cmd="ls {}".format(folder)
        return sp.getoutput(cmd).split("\n")

def get_basecall_folder(wildcards):
    folder="/".join([INDIR,basecalled_bam_string.format(wildcards=wildcards)])
    return folder

def get_clair_model(wildcards):
    my_flowcell = samples.loc[wildcards.SAMPLEID,"Flowcell"]
    if my_flowcell=="R9":
        return '{}/clair3_models/r941_prom_sup_g5014'.format(config["clairmodelpath"])
    return '{}/rerio/clair3_models/r1041_e82_400bps_sup_v500'.format(config["clairmodelpath"])


def get_hpdp_png_names(wildcards):
    if wildcards.STRATEGY=="RU":
        bedfile="".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    else:
        bedfile="workflow/resources/hpdp_targets.bed"
    df=pd.read_table(bedfile, sep="\t", header=None)
    regionnames=df[0].tolist()
    return ["".join([PREFIX_REGEX, ".hp_dp.detail_plot.png.", x, ".png"]) for x in regionnames]

def get_target_genes(wildcards):
    genes=samples.loc[wildcards.SAMPLEID, "TargetGene"].split("_")
    return genes

def get_sv_outputs(wildcards):
    callers=["None"]
    if config["all_targets"]:
        callers+=["cuteSV", "svim", "sniffles"]
        return callers
    elif config["outputs"]["cuteSV"]:
        callers+=["cutesv"]
    if config["outputs"]["svim"]:
        callers+=["svim"]
    if config["outputs"]["sniffles"]:
        callers+=["sniffles"]
    return callers

def get_report_inputs(wildcards):
    endings=[]
    if wildcards.STRATEGY=="RU":
        endings+=["target.plot_readlengths.png"]
        if config["outputs"]["VEP"] or config["allTargets"]:
            endings+=["vep.target.snippet.html"]
        if config["outputs"]["clair3"] or config["allTargets"]:
            endings+=["target.plot_indel_quality.png", "target.plot_snv_quality.png"]
    if config["outputs"]["alignBam"] or config["allTargets"]:
        endings+=["plot_readlengths.png", "plot_depth_coverage.png", "read_efficiency.png", "yield_efficiency.png", "n50_efficiency.png"]
    if config["outputs"]["clair3"] or config["allTargets"]:
        endings+=["plot_indel_quality.png", "plot_snv_quality.png", "clair3.notPhased.snpPlot.png"]
    if config["outputs"]["CNVcalls"] or config["allTargets"]:
       endings+=["called_cnv.detail_plot.png"]
    if config["outputs"]["VEP"] or config["allTargets"]:
        endings+=["vep.pathogenic.snippet.html"]
    if config["outputs"]["phaseQC"] or config["allTargets"]:
        endings+=["clair3.phased.whatshap_plot.png", "hp_dp_long_complete.txt", "hp_dp.snippet.html"]
    return [ ".".join([PREFIX_REGEX, x]) for x in endings]

def get_target_bed(wildcards):
    if wildcards.STRATEGY == "RU":
        return "".join([config["bedfiledir"],"/",samples.loc[wildcards.SAMPLEID, "BedFile"]])
    else:
        return "workflow/resources/hpdp_targets.bed"

