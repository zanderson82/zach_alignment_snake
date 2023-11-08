import pandas as pd
import subprocess as sp
import numpy as np


def sort_snps(x):
    if x <1:
        return "SNP"
    return "INDEL"

#run a few shell commands to filter vep with its own built in tool, get the length of the header, and pull the names of additional columns in the INFO field
vep_filter_cmd=" ".join(["filter_vep -i", snakemake.input.vep_vcf, "-o", snakemake.output.vep_intermediate, """--filter "(MAX_AF <= 0.01 or not MAX_AF) and (gnomADg_AF <= 0.01 or not gnomADg_AF) and (FILTER != LowQual)" --format vcf --force_overwrite"""])
vepFilterOut=sp.getoutput(vep_filter_cmd)
skiplines= sp.getoutput(" ".join(["grep '^##*'", snakemake.output.vep_intermediate,"| wc -l"]))
skiplines=int(skiplines)
info_header=sp.getoutput(" ".join(["grep '^##INFO=<ID=CSQ'", snakemake.output.vep_intermediate, "| cut -d ' ' -f 7"]))
info_header = info_header.split('">')[0]

# read filtered vep file in
df = pd.read_table(snakemake.output.vep_intermediate, sep="\t", skiprows=skiplines, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
df["HP"] = df.SAMPLE.apply(lambda x: x.split(":")[0])
df["DP_ALL"] = df.SAMPLE.apply(lambda x: int(x.split(":")[2]))
df["IGV"] = df.CHROM.str.cat(df.POS.astype(str), sep=":")

df_phased = df[df.SAMPLE.apply(lambda x: len(x.split(":"))==6)]
df_unphased = df[df.SAMPLE.apply(lambda x: len(x.split(":"))!=6)]


df_phased["AF"] = df_phased.SAMPLE.apply(lambda x: float(x.split(":")[4]))
df_phased["DP_HP1"] = df_phased.DP_ALL.astype("float") * df_phased.AF
df_phased.DP_HP1 = df_phased.DP_HP1.apply(np.round)
df_phased.DP_HP1 = df_phased.DP_HP1.astype(int)
df_phased["DP_HP2"] = df_phased.DP_ALL - df_phased.DP_HP1
df_phased["PHASEBLOCK"] = df_phased.SAMPLE.apply(lambda x: x.split(":")[-1])

snp_lengths = df_phased.REF.apply(len) - df_phased.ALT.apply(len)
df_phased["TYPE"] = snp_lengths.apply(sort_snps)
call_sources = df_phased.INFO.apply(lambda x: x.split(";")[0])
df_phased["SOURCE"] = call_sources.str.replace("F", "Full Alignment").str.replace("P", "Pileup")
info_field = df_phased.INFO.apply(lambda x: x.split("CSQ=")[-1])
df_info = info_field.str.split("|", expand=True)
df_info.columns = info_header.split("|")
df_out_phased = df_phased[["CHROM", "POS", "IGV", "REF", "ALT", "QUAL", "TYPE", "HP", "DP_ALL", "DP_HP1", "DP_HP2", "PHASEBLOCK", "SOURCE"]].join(df_info, how="inner", lsuffix="CALL")
df_out_phased.to_csv(snakemake.output.vep_lt1_phased_vcf, index=False)

df_unphased["ALT1"] = df_unphased.ALT.apply(lambda x: x.split(",")[0])
df_unphased["ALT2"] = df_unphased.ALT.apply(lambda x: x.split(",")[-1])
df_unphased.loc[df_unphased.ALT1==df_unphased.ALT2,"ALT2"] = "-"
df_unphased["AF_1"] = df_unphased["SAMPLE"].apply(lambda x: float(x.split(":")[-1].split(",")[0]))
df_unphased["AF_2"] = df_unphased["SAMPLE"].apply(lambda x: float(x.split(":")[-1].split(",")[-1]))
df_unphased.loc[df_unphased.AF_1==df_unphased.AF_2,"AF_2"] = 0

df_unphased["DP_HP1"] = df_unphased.DP_ALL.astype("float") * df_unphased.AF_1
df_unphased.DP_HP1 = df_unphased.DP_HP1.apply(np.round)
df_unphased.DP_HP1 = df_unphased.DP_HP1.astype(int)
df_unphased["DP_HP2"] = df_unphased.DP_ALL.astype("float") * df_unphased.AF_2
df_unphased.DP_HP2 = df_unphased.DP_HP2.apply(np.round)
df_unphased.DP_HP2 = df_unphased.DP_HP2.astype(int)
df_unphased["DP_HP3"] = df_unphased.DP_ALL - df_unphased.DP_HP1 - df_unphased.DP_HP2

snp_lengths1 = df_unphased.REF.apply(len) - df_unphased.ALT1.apply(len)
df_unphased["TYPE1"] = snp_lengths1.apply(sort_snps)
snp_lengths2 = df_unphased.REF.apply(len) - df_unphased.ALT2.apply(len)
df_unphased["TYPE2"] = snp_lengths2.apply(sort_snps)
df_unphased.loc[df_unphased.TYPE1==df_unphased.TYPE2,"TYPE2"] = "-"
call_sources = df_unphased.INFO.apply(lambda x: x.split(";")[0])
df_unphased["SOURCE"] = call_sources.str.replace("F", "Full Alignment").str.replace("P", "Pileup")
info_field = df_unphased.INFO.apply(lambda x: x.split("CSQ=")[-1])
df_info = info_field.str.split("|", expand=True)
df_info.columns = info_header.split("|")
df_out_unphased = df_unphased[["CHROM", "POS", "IGV", "REF", "ALT1", "ALT2", "QUAL", "TYPE1", "TYPE2", "HP", "DP_ALL", "DP_HP1", "DP_HP2", "DP_HP3", "SOURCE"]].join(df_info, how="inner", lsuffix="CALL")
df_out_unphased.to_csv(snakemake.output.vep_lt1_notPhased_vcf, index=False)