#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

invlength <- as.numeric(args[1])
setwd(paste0("~/Dropbox/projects/branch-durations/slim-simulations/", args[2]))

suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(invClust))
suppressPackageStartupMessages(library(tidyverse))

vcf <- readVcf(paste0("samples.vcf"))
mat <- genotypeToSnpMatrix(vcf)
geno <- mat$genotypes
annot <- read.table("relate.mut", se=";", header=TRUE) 
annot <- annot %>%
  mutate("chromosome" = 0) %>%
  rename(`snp.name` = snp) %>%
  rename(position = pos_of_snp) %>%
  select(chromosome, `snp.name`, position)

truth <- read.table("carriers.txt", header=FALSE, sep=",")

# roi1<-data.frame(chr=0, LBP=2500000 - as.integer(invlength/2), RBP=2500000 + as.integer(invlength/2), reg= "inv1")
# roi2<-data.frame(chr=0, LBP=2500000 - as.integer(invlength), RBP=2500000 + as.integer(invlength), reg= "inv2")
# roi3<-data.frame(chr=0, LBP=1000000 - as.integer(invlength/2), RBP=1000000 + as.integer(invlength/2), reg= "inv3")
# roi4<-data.frame(chr=0, LBP=4000000 - as.integer(invlength/2), RBP=4000000 + as.integer(invlength/2), reg= "inv4")

left = c(2500000 - as.integer(invlength/2),
         2500000 - as.integer(invlength),
         1000000 - as.integer(invlength/2),
         4000000 - as.integer(invlength/2))
right = c(2500000 + as.integer(invlength/2),
          2500000 + as.integer(invlength),
          1000000 + as.integer(invlength/2),
          4000000 + as.integer(invlength/2))

results <- data.frame(inferred_frequency = c(), inferred_genotypes = c())
for(i in 1:length(left)) {
  roi <- data.frame(chr=0, LBP=left[i], RBP=right[i], reg= "inv1")
  invcall1<-invClust(roi=roi, wh = 1, geno=geno, annot=annot, dim=2)
  invcall1_pred <- data.frame(invGenotypes(invcall1)) %>%
    rename(t = invGenotypes.invcall1.) %>%
    mutate(total = ifelse(t == "NI/NI", 0, ifelse(t == "NI/I", 0.5, 1)))
  check1 <- sum(invcall1_pred$total)
  check2 <- cor(invcall1_pred$total, truth$V2)^2
  results <- rbind(results, data.frame(inferred_frequency = c(check1), inferred_genotypes = c(check2)))
}

write.table(results, file = "invclust.txt", quote=FALSE, row.names=FALSE)
