#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]

outputparts=stringr::str_split(output, "\\.")[[1]]
mode=outputparts[length(outputparts)]

suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(kableExtra)))


df <- read.csv(input, header=TRUE)

df %>% subset(select=-c(File)) %>% kbl() %>%
    column_spec(3, color="white", background=spec_color(df$Depth, scale_from=c(0,45))) %>%
    column_spec(6, color="white", background=spec_color(df$PercentAssigned, scale_from=c(0,100))) %>%
    column_spec(7, color=spec_color(df$Phaseblocks, scale_from=c(0,10), direction=1, option="A")) %>% save_kable(file=output, self_contained=T)

df %>% subset(select=-c(File)) %>% kbl() %>% save_kable(file=output, self_contained=T)