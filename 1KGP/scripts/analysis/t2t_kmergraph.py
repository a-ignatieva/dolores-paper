#!/usr/bin/env python
# coding: utf-8

import sys
import os
import gzip
import glob
import math
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
import tqdm
import matplotlib.patches as patches
from matplotlib.lines import Line2D

sys.path.append("..")

ann_dir = "t2t/annotations/"
kmer_dir = "t2t/kmers/"
output_dir = "t2t/plots/"

band = sys.argv[1]
startgene = sys.argv[2]
outname = str(band) + "_kmergraph.pdf"

dna_dict = {"A": "T", "C": "G", "T": "A", "G": "C"}

def reverse_seq(string):
    new_string = ""
    for s in string:
        new_string += dna_dict[s]
    return new_string

def flip_seq(string):
    new_string = string[::-1]
    return new_string

p1 = p2 = s = None
with open("files/liftover.csv", "r") as file:
    for line in file:
        line = line.strip().split(",")
        if line[0] == band:
            p1 = int(line[2])
            p2 = int(line[3])
            break
    if p1 is not None:
        with gzip.open(ann_dir + "Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz", "rt") as gfile:
            for l, line in enumerate(gfile):
                line = line.strip()
                if l > 4:
                    line = line.split("\t")
                    if line[2] == "gene":
                        s = int(line[3])
                        details_ = line[-1].split(";")
                        details = {}
                        for d in details_:
                            d = [dd.strip() for dd in d.split('"')]
                            if len(d) > 1 and d[0] not in details:
                                details[d[0]] = d[1]
                        if "gene_name" in details:
                            gene_name = details["gene_name"]
                            if gene_name == startgene:
                                p1 = p1 - s
                                p2 = p2 - s
                                break
print(startgene, s, p1, p2)

with open("t2t/plots/counts.csv", "r") as file:
    i = 0
    IDs = {}
    carriers = {}
    for l, line in enumerate(file):
        line = line.strip().split(",")
        if l == 1:
            for ll in line:
                if ll != "":
                    IDs[ll] = i
                    i += 1
        elif l > 2:
            carriers[line[0]] = [float(ll) for ll in line[1:] if ll != ""]

fig, axs = plt.subplots(
    1,
    1,
    figsize=(12, 20),
    gridspec_kw={
        "height_ratios": [1.0]
    },
)

directory = os.fsencode(output_dir)
filelist = glob.glob(os.path.join(kmer_dir, band + '_kmer_counts_*.csv.gz'))
kmers_ref = set()
kmers_ref_list = {}

L = 0
M = 1e10
with gzip.open(kmer_dir + band + "_kmer_counts_T2T-CHM13v2.none.csv.gz", "rt") as file:
    for line in file:
        line = line.strip().split(",")
        kmers_ref.add(line[0])
        pos = int(line[2].split(";")[0])
        L = max(L, pos)
total_kmers = len(kmers_ref)
print("Total kmers:", total_kmers)
prevpos = 0
if p1 is None:
    d = int(L/50)
    with gzip.open(kmer_dir + band + "_kmer_counts_T2T-CHM13v2.none.csv.gz", "rt") as file:
        for line in file:
            line = line.strip().split(",")
            kmer = line[0]
            kmer_ = flip_seq(reverse_seq(kmer))
            pos = int(line[2].split(";")[0])
            if int(line[1]) == 1 and kmer_ not in kmers_ref and abs(pos - prevpos) >= d:
                kmers_ref_list[kmer] = pos
                prevpos = pos
                M = min(M, pos)
else:
    L = 0
    d = (p2-p1)*3/50
    with gzip.open(kmer_dir + band + "_kmer_counts_T2T-CHM13v2.none.csv.gz", "rt") as file:
        for line in file:
            line = line.strip().split(",")
            kmer = line[0]
            kmer_ = flip_seq(reverse_seq(kmer))
            pos = int(line[2].split(";")[0])
            if p1 - (p2 - p1) <= pos <= p2 + (p2 - p1):
                if int(line[1]) == 1 and kmer_ not in kmers_ref and abs(pos - prevpos) >= d:
                    kmers_ref_list[kmer] = pos
                    prevpos = pos
                    M = min(M, pos)
                    L = max(L, pos)
kmers_ref = set(kmer for kmer in kmers_ref_list)
unique_kmers = len(kmers_ref)
w = (L-M) / len(kmers_ref)
print("Total unique kmers:", unique_kmers)

axs.hlines(
    xmin=p1,
    xmax=p2,
    y=1,
    linewidth=5,
    color="black",
)
yy = 0
inds = [v for v in kmers_ref_list.values()]
kmers_to_label = [x for _, x in sorted(zip(inds, [k for k in kmers_ref_list]))]
clrs = sns.color_palette("husl", n_colors=len(kmers_to_label))
col_pal = {kmers_to_label[i]: clrs[i] for i in range(len(kmers_to_label))}
with tqdm.tqdm(total = len(filelist)) as pbar:
    for i, filename in enumerate(reversed(sorted(filelist))):
        kmers_seq = set()
        filename_ = filename.strip().split("_")[3].split(".")
        name = filename_[0]
        mp = filename_[1]
        with gzip.open(filename, "rt") as file:
            n = name + "." + mp
            if n in IDs:
                reflabel_ = n + " " + str(carriers[band][IDs[n]]) + " "
            else:
                reflabel_ = n + " -- "
            axs.text(
                M,
                yy,
                reflabel_,
                ha="right",
                va="center",
                fontsize=10,
            )
            for line in file:
                line = line.strip().split(",")
                kmer = line[0]
                kmer_ = flip_seq(reverse_seq(kmer))
                ct = int(line[1])
                p = [int(pp) for pp in line[2].split(";")]
                if kmer in kmers_ref:
                    for pos in p:
                        # axs.scatter(
                        #     x = pos,
                        #     y = yy,
                        #     color = col_pal[kmer],
                        #     s = 40,
                        #     marker="s",
                        # )
                        axs.hlines(
                            xmin = pos - w/2,
                            xmax = pos + w/2,
                            y = yy,
                            linewidth=10,
                            color = col_pal[kmer],
                        )
                        axs.text(
                            pos,
                            yy,
                            ">",
                            ha = "center",
                            va = "center",
                            fontsize = 10,
                        )
                if kmer_ in kmers_ref:
                    for pos in p:
                        axs.hlines(
                            xmin = pos - w/2,
                            xmax = pos + w/2,
                            y = yy,
                            linewidth=10,
                            color = col_pal[kmer_],
                        )
                        axs.text(
                            pos,
                            yy,
                            "<",
                            ha="center",
                            va="center",
                            fontsize=10,
                        )
        yy -= 1
        pbar.update(1)

plt.axis("off")
plt.tight_layout()
plt.savefig(output_dir + "/" + outname, dpi=300, bbox_inches="tight")




