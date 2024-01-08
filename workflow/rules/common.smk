import pandas as pd
import subprocess as sp
import glob

samples = pd.read_table(config["samples"], sep="\t").set_index('SampleID')

THREADS=config["threads"]
#FLOWCELL=config["flowcell"]
REFGENOME=config["refgenome"]
ONTMMIFILE=config["ontmmifile"]
WORKDIR=config["working_dir"]
FINALDIR=config["final_dir"]
INDIR=config["input_dir"]
PROJECT=config["project"]
#STRATEGY=config["strategy"]
BEDFILEDIR=config["bedfiledir"]

PREFIX_REGEX="{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}_{GENE}-{METH}-{MB}/{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}_{GENE}-{METH}-{MB}"
SAMPLE_WORKPATH="".join([WORKDIR, "/", PREFIX_REGEX]) # SAMPLE_WORKPATH should be used in all inputs up to moving the files to /n/alignments
LOG_REGEX="{SAMPLEID}-NP-{STRATEGY}-{PROJECT_ID}-{OUTSIDE_ID}_{GENE}-{METH}-{MB}"

def get_flowcell(wildcards):
    sampleID=wildcards.SAMPLEID
    return samples.loc[sampleID,"Flowcell"]

def apply_suffix(wildcards, ending, sampleID):
    # applies a suffix to the path of a sample folder and sample prefix, e.g. M1079-NP-{STRATEGY}-Project-OutsideID-pb/M1079-NP-{STRATEGY}-Project-OutsideID-pb.vcf.gz
    return ".".join([get_output_name(wildcards, sampleID), ending])

def get_output_name(wildcards, sampleID):
    #returns the path to a sample folder and sample prefix, e.g. M1079-NP-{STRATEGY}-Project-OutsideID-pb/M1079-NP-{STRATEGY}-ProjectOutsideID-pb
    return "{sid}-NP-{strat}-{project}-{oid}_{gene}-{meth}-{member}/{sid}-NP-{strat}-{project}-{oid}_{gene}-{meth}-{member}".format(sid=sampleID, strat=samples.loc[sampleID, "Strategy"], project=config["project"], oid=samples.loc[sampleID,"ExternalID"], gene=samples.loc[sampleID,"TargetGene"], meth=samples.loc[sampleID,"Methylation"],member=samples.loc[sampleID,"Member"])

def get_output_dir(wildcards):
    # returns the working (Franklin) directory with specific folder, an absolute path. e.g. /data/alignments/M1079-NP-{STRATEGY}-Project-OutsideID-pb/
    sampleID=wildcards.SAMPLEID
    return "{outdir}/{sid}-NP-{strat}-{project}-{oid}_{gene}-{meth}-{member}".format(outdir=WORKDIR, sid=sampleID, strat=samples.loc[sampleID, "Strategy"], project=config["project"], oid=samples.loc[sampleID,"ExternalID"], gene=samples.loc[sampleID,"TargetGene"], meth=samples.loc[sampleID,"Methylation"], member=samples.loc[sampleID,"Member"])

def get_final_dir(wildcards):
    # returns the destination (McClintock) directory with specific folder, a relative path
    sampleID=wildcards.SAMPLEID
    return "{finaldir}/{sid}-NP-{strat}-{project}-{oid}_{gene}-{meth}-{member}".format(finaldir=FINALDIR, sid=sampleID, strat=samples.loc[sampleID, "Strategy"], project=config["project"], oid=samples.loc[sampleID,"ExternalID"], gene=samples.loc[sampleID,"TargetGene"], meth=samples.loc[sampleID,"Methylation"], member=samples.loc[sampleID,"Member"])

def get_only_multisample_targets(wildcards):
    f = open(config["targetfile"], "r")
    targetsamples = f.read().split("\n")
    f.close()
    final_targets = []
    for ts in targetsamples:
        all_targets = get_trio_files(wildcards, ts, targetsamples)
        family_dir=get_final_dir(wildcards, ts)
        final_targets += ["/".join([family_dir,x]) for x in all_targets]
    return final_targets


def get_final_targets(wildcards):
    f = open(config["targetfile"], "r")
    targetsamples = f.read().split("\n")
    f.close()
    final_targets = []
    for ts in targetsamples:
        strategy=samples.loc[ts,"Strategy"]
        file_endings= ["phased.bam", "phased.bam.bai"] #phased bam
        file_endings += ["clair3.phased.vcf.gz", "clair3.notPhased.vcf.gz", "clair3.phased.vcf.gz.tbi", "clair3.notPhased.vcf.gz.tbi"] #clair stuff
        file_endings += ["sv_cutesv.notPhased.vcf", "sv_sniffles.notPhased.vcf", "sv_svim.notPhased.vcf"] #svs
        file_endings += ["sv_cutesv.phased.vcf", "sv_sniffles.phased.vcf", "sv_svim.phased.vcf"] # phased svs
        file_endings += ["clair3.phased.vep.vcf", "clair3.phased.vep.af_lt_1.csv"] #vep SeqFirst project not using VEP
        if strategy == "RU":
            file_endings += ["clair3.phased.phasing_stats.tsv", "phased.target.cramino.stats", "phased.target.bam", "phased.target.bam.bai", "phased.target.coverage.tsv"]
        else:    
            file_endings += ["clair3.phased.phasing_stats.tsv", "phased.cramino.stats"] # alignment QC stats
        all_targets = [apply_suffix(wildcards, x, ts) for x in file_endings] #add trio stuff
        #all_targets += get_trio_files(wildcards, ts, targetsamples) # this checks if trio files should be made, and if so adds them in the subfolder FAMILY_multisample
        #family_dir=get_final_dir(wildcards, ts) # add the family folder to the front of each file.
        final_targets += all_targets
    return final_targets

# this version targets all samples in the reference TSV.
def get_final_targets_all(wildcards):
    allSamples=samples.index.tolist()
    final_targets = []
    for ts in allSamples:
        strategy=samples.loc[ts,"Strategy"]
        file_endings= ["phased.bam", "phased.bam.bai"] #phased bam
        file_endings += ["clair3.phased.vcf.gz", "clair3.notPhased.vcf.gz", "clair3.phased.vcf.gz.tbi", "clair3.notPhased.vcf.gz.tbi"] #clair stuff
        file_endings += ["sv_cutesv.notPhased.vcf", "sv_sniffles.notPhased.vcf", "sv_svim.notPhased.vcf"] #svs
        file_endings += ["sv_cutesv.phased.vcf", "sv_sniffles.phased.vcf", "sv_svim.phased.vcf"] # phased svs
        file_endings += ["clair3.phased.vep.vcf", "clair3.phased.vep.af_lt_1.csv"] #vep SeqFirst project not using VEP
        if strategy == "RU":
            file_endings += ["clair3.phased.phasing_stats.tsv", "phased.target.cramino.stats", "phased.target.bam", "phased.target.bam.bai", "phased.target.coverage.tsv"]
        else:    
            file_endings += ["clair3.phased.phasing_stats.tsv", "phased.cramino.stats"]
        all_targets = [apply_suffix(wildcards, x, ts) for x in file_endings] #add trio stuff
        #all_targets += get_trio_files(wildcards, ts, targetsamples) # this checks if trio files should be made, and if so adds them in the subfolder FAMILY_multisample
        #family_dir=get_final_dir(wildcards, ts) # add the family folder to the front of each file.
        final_targets += all_targets
    return final_targets

def get_target_bams(wildcards):
    #return config["libraries"][wildcards.SAMPLEID].split(" ")
    cmd="".join(["ls ", INDIR, "/", wildcards.SAMPLEID, "*"])
    return sp.getoutput(cmd).split("\n")

def get_clair_model(wildcards):
    my_flowcell = get_flowcell(wildcards)
    if my_flowcell=="R9":
        return '{}/clair3_models/r941_prom_sup_g5014'.format(config["clairmodelpath"])
    return '{}/rerio/clair3_models/r1041_e82_400bps_sup_v420'.format(config["clairmodelpath"])
