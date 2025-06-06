#!/usr/bin/env python
# coding: utf-8

import sys
import os

import vcf
import gzip

import numpy as np
from scipy.stats import pearsonr
import math

from liftover import get_lifter

converter = get_lifter("hg19", "hg38")

sys.path.append("..")

results_file = "files/clades_pvalues_all_events_annotated_varNe.csv"
nodes_file = "files/bands_to_freq_nodes.csv"
genotypes_file = "1KGP_data/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz"
samples_file = "files/HPRC_samples.csv"
output_file = "files/counts_HPRC_updated.csv"

samples = {}
with open(samples_file, "r") as file:
    print("Reading in samples file...")
    for line in file:
        line = line.strip().split(",")
        samples[line[0]] = [0, 0, 0, 0, line[1]]
print(len(samples))

bands = {}
with open(results_file, "r") as file:
    print("Reading in results file...")
    next(file)
    for line in file:
        line = line.strip().split(sep=",")
        bands[line[21]] = int(line[0][3:])

with open(output_file, "w") as ofile:
    ofile.write(",")
    for sample in samples:
        ofile.write(sample + ".pat," + sample + ".mat,")
    ofile.write("\n,")
    for sample in samples:
        ofile.write(samples[sample][4] + "," + samples[sample][4] + ",")
    ofile.write("\n")

    for band, chromosome in bands.items():

        freq_file = "results_1000GP/chr" + band + "_freq.csv"
        relate_file = (
            "trees/relate_1000GP/ancmut/1000GP_Phase3_mask_prene_chr"
            + str(chromosome)
            + ".mut.gz"
        )

        for sample in samples:
            samples[sample] = [0, 0, 0, 0, 0]

        nodes = {}
        with open(nodes_file, "r") as file:
            print("Reading in nodes file...")
            for line in file:
                line = line.strip().split(",")
                nodes[line[0]] = int(line[1])
        print(len(nodes))

        with open(results_file, "r") as file:
            print("Reading in results file...")
            for line in file:
                line = line.strip().split(sep=",")
                if line[21] == band:
                    pop = line[1]
                    start = converter[str(chromosome)][int(line[10])][0][1]
                    end = converter[str(chromosome)][int(line[11])][0][1]

        SNPs_hg19 = set()
        with open(freq_file, "r") as file:
            print("Reading in frequencies file...")
            for line in file:
                line = line.strip().split(sep=",")
                if int(line[0]) == nodes[band]:
                    SNPs_hg19.add(int(line[12]))

        SNPs = {}
        with gzip.open(relate_file, "rt") as file:
            print("Reading in relate file...")
            next(file)
            for line in file:
                line = line.strip().split(sep=";")
                if int(line[1]) in SNPs_hg19:
                    p = converter[str(chromosome)][int(line[1])][0]
                    print(int(line[1]), p, int(p[1]))
                    SNPs[int(p[1])] = (line[10][0], line[10][2], int(line[1]), line[3])  # hg38 pos: (ancestral, alt, hg19 pos, rsid)
        print(SNPs)

        vcf_reader = vcf.Reader(
            filename=genotypes_file,
            compressed=True,
        )
        records = vcf_reader.fetch("chr" + str(chromosome), start, end)

        samples_found = set()
        snp_count = 0
        for record in records:
            if record.POS in SNPs:
                print(record, SNPs[record.POS])
                sample_count = 0
                snp_count += 1
                if record.REF not in SNPs[record.POS]:
                    print("REF do not match")
                    continue
                if record.REF == SNPs[record.POS][0]:
                    carrier = 1
                else:
                    carrier = 0
                # print("carrier:", carrier)
                for sample in record.samples:
                    if sample.sample in samples:
                        samples_found.add(sample.sample)
                        sample_count += 1
                        if sample["GT"][0] != ".":
                            samples[sample.sample][0] += int(int(sample["GT"][0]) == carrier)
                            samples[sample.sample][1] += 1
                        if sample["GT"][2] != ".":
                            samples[sample.sample][2] += int(int(sample["GT"][2]) == carrier)
                            samples[sample.sample][3] += 1
                        # print(int(int(sample["GT"][0]) == carrier), int(int(sample["GT"][2]) == carrier))
                print("samples found:", sample_count)

        print("SNPs found:", snp_count, "out of", len(SNPs))
        for sample in samples:
            if samples[sample][1] == 0:
                samples[sample][0] = "NA"
            else:
                samples[sample][0] = int(samples[sample][0]/samples[sample][1] > 0.5)
            if samples[sample][3] == 0:
                samples[sample][2] = "NA"
            else:
                samples[sample][2] = int(samples[sample][2]/samples[sample][3] > 0.5)

        outstring = band
        for sample, (pat, _, mat, _, _) in samples.items():
            outstring += "," + str(pat) + "," + str(mat)
        outstring += "\n"
        ofile.write(outstring)
