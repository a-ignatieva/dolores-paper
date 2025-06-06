#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pickle

import tskit
import tszip
import stdpopsim
import vcf
import pysam

import numpy as np
from scipy.stats import pearsonr
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

sys.path.append("..")
import argbd.cladedurations_
import argbd.simulations
import argbd.viz

col_green = "#228833"
col_red = "#EE6677"
col_purp = "#AA3377"
col_blue = "#66CCEE"
col_yellow = "#CCBB44"
col_indigo = "#4477AA"
col_grey = "#BBBBBB"


chromosome = sys.argv[1]
band = sys.argv[2]
name = chromosome + band + "_time"
clades_to_time = sys.argv[3]
clades_to_include = sys.argv[4]
plot_mrca = True
time_all = False
if clades_to_time != "None":
    if clades_to_time == "-1":
        clades_to_time = []
        time_all = True
    else:
        clades_to_time = [int(i) for i in clades_to_time.split(sep=",")]
else:
    clades_to_time = []
if clades_to_include != "None":
    clades_to_include = [int(i) for i in clades_to_include.split(sep=",")]
else:
    clades_to_include = []

def get_bitset(n_set, samp_to_ind):
    e = [0] * int(len(samp_to_ind) / 2)
    for n in n_set:
        e[samp_to_ind[n]] += 1
    return e


output_dir = "results_1000GP"
mask_loc = "files/20140520.pilot_mask.autosomes.bed"
results_file = (
    "files/clades_pvalues_all_events_annotated.csv"
)
superdups_file = "files/genomicSuperDups.txt"
svs_file = "files/common_1000g.bed"
annot_file = "files/ncbiRefSeqCurated.txt"
svgenotypes_file = "files/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"
trees_loc = "trees/relate_1000GP/trees/chunks/"
time_file = "results_1000GP/" + chromosome + band + "_time.csv"
mutmap_file = "results_1000GP/" + chromosome + band + "_mutmap.csv"

mask = argbd.simulations.read_genome_mask_bed(mask_loc, chromosome)
(
    chromosome_lengths,
    chromosome_starts,
    centromeres,
    centromeres_mask,
) = argbd.simulations.read_genome_info(output_dir + "/chromosomes_info.txt")
colors = {
    "EUR": col_indigo,
    "AFR": "orange",
    "EAS": col_green,
    "SAS": col_purp,
    "AMR": col_blue,
    "ALL": "black"
}
num_samples = {
    "EUR": 1006,
    "AFR": 1322,
    "EAS": 1008,
    "SAS": 978,
    "AMR": 694,
    "ASIA": 1008,
    "ADMIX": 694,
}

species = stdpopsim.get_species("HomSap")
contig = species.get_contig(chromosome=chromosome, genetic_map="HapMapII_GRCh37")
recombination_map = contig.recombination_map

md, pops, groups = argbd.simulations.read_in_pop_metadata(
    "trees/relate_1000GP/data/1000GP_Phase3.poplabels"
)
samp_to_ind_all = {}
inds = {}
i = 0
for key, val in md.items():
    if val["ID"] not in inds:
        inds[val["ID"]] = i
        i += 1
    samp_to_ind_all[key] = inds[val["ID"]]

fig, ax_time = plt.subplots(
    1,
    1,
    figsize=(8, 3),
    gridspec_kw={"height_ratios": [1.0]},
)

xmin = math.inf
xmax = -1
ymax = 0
bitset_all = [0] * int(len(samp_to_ind_all) / 2)
chunkmin = math.inf

if os.path.exists(time_file):
    with open(results_file, "r") as file:
        print("Reading in time file...")
        if not time_all:
            for line in file:
                line = line.strip().split(sep=",")
                if line[0] == chromosome:
                    if line[21] == chromosome[3:] + band:
                        node_id = int(line[20])
                        if node_id in clades_to_time:
                            pop = line[1]
                            pop_ = "_" + pop
                            start = int(line[10])
                            end = int(line[11])
                            clade_id = int(line[5])
                            left_mut = int(line[14])
                            right_mut = int(line[15])
                            chunkindex = int(line[18])
                            treeindex = int(line[19])
                            node_id = int(line[20])
                            cladesize = int(line[8])
                            p1 = float(line[6])
                            p2 = float(line[7])
                            chunkmin = min(chunkmin, chunkindex)
                            samp_to_ind = pickle.load(
                                open(
                                    trees_loc
                                    + chromosome
                                    + "/1000GP_Phase3_mask_prene_"
                                    + chromosome
                                    + pop_
                                    + "_samp_to_ind.pickle",
                                    "rb",
                                )
                            )
                            ts = tszip.decompress(
                                trees_loc
                                + chromosome
                                + "/1000GP_Phase3_mask_prene_"
                                + chromosome
                                + "_chunk"
                                + str(chunkindex)
                                + pop_
                                + ".trees.tsz"
                            )
                            mapp = pickle.load(
                                open(
                                    trees_loc
                                    + chromosome
                                    + "/1000GP_Phase3_mask_prene_"
                                    + chromosome
                                    + pop_
                                    + "_map.pickle",
                                    "rb",
                                )
                            )

                            t = ts.at_index(treeindex)
                            samples = {s for s in t.samples(node_id)}
                            bitset_true = get_bitset(samples, samp_to_ind)
                            if node_id in clades_to_include:
                                for s in samples:
                                    ss = np.where(mapp == s)[0][0]
                                    bitset_all[samp_to_ind_all[ss]] += 1

                            xmin = min(xmin, start)
                            xmax = max(xmax, end)

                            with open(time_file, "r") as tfile:
                                for tline in tfile:
                                    tline = tline.strip().split(sep=",")
                                    node_id_ = int(tline[1])
                                    if node_id == node_id_:
                                        r2 = float(tline[5])
                                        left = float(tline[19])
                                        right = float(tline[20])
                                        t1 = float(tline[21])
                                        t2 = float(tline[22])
                                        if r2 > 0.90:
                                            if r2 >= 0.95:
                                                c = colors[pop]
                                            else:
                                                c = "grey"
                                            p = patches.Rectangle(
                                                (left, t1),
                                                right - left,
                                                t2-t1,
                                                facecolor=c,
                                            )
                                            ymax = max(ymax, t2)
                                            ax_time.add_patch(p)
        else:
            pop = "ALL"
            with open(time_file, "r") as tfile:
                for tline in tfile:
                    tline = tline.split(sep=",")
                    if tline[0] == "ALL":
                        r2 = float(tline[5])
                        left = float(tline[19])
                        right = float(tline[20])
                        t1 = float(tline[21])
                        t2 = float(tline[22])
                        if r2 > 0.90:
                            if r2 >= 0.95:
                                c = colors[pop]
                            else:
                                c = "grey"
                            p = patches.Rectangle(
                                (left, t1),
                                right - left,
                                t2 - t1,
                                facecolor=c,
                            )
                            ymax = max(ymax, t2)
                            ax_time.add_patch(p)

if plot_mrca:
    print("Getting ts at chunk", chunkmin, "at position", xmin)
    ts = tszip.decompress(
        trees_loc
        + chromosome
        + "/1000GP_Phase3_mask_prene_"
        + chromosome
        + "_chunk"
        + str(chunkmin)
        + ".trees.tsz"
    )
    t = ts.at(xmin)
    while t.interval[1] < xmax:
        if t.index > 0 and t.index < ts.num_trees - 1:
            ax_time.hlines(xmin=t.interval[0], xmax=t.interval[1], y=t.time(t.root), color="black")
        if t.index == ts.num_trees - 2:
            chunkmin += 1
            ts = tszip.decompress(
                trees_loc
                + chromosome
                + "/1000GP_Phase3_mask_prene_"
                + chromosome
                + "_chunk"
                + str(chunkmin)
                + ".trees.tsz"
            )
            t = ts.first()
        t.next()


xmin = xmin - 0.02 * (xmax - xmin)
xmax = xmax + 0.02 * (xmax - xmin)

xticks = [i*100000 for i in range(0, 3000) if xmin <= i*100000 <= xmax]
xlabels = [str(i/1000000) for i in xticks]
plt.setp(ax_time, xticks=xticks, xticklabels=xlabels)
xticks = [i*10000 for i in range(0, 30000) if xmin <= i*10000 <= xmax]
ax_time.scatter(x=xticks, y = [0]*len(xticks), marker="|", zorder=10, clip_on=False, color="black", s=2)

ax_time.set_title("Age estimate")
ax_time.set_xlim([xmin, xmax])
ax_time.set_ylim([0, ymax*1.05])
ax_time.set_ylabel("Time (generations)")
ax_time.set_xlabel("Genome position (Mb) on " + chromosome)

plt.tight_layout()
plt.savefig(output_dir + "/" + name + ".png", dpi=300, bbox_inches="tight", pad_inches=0.1)
