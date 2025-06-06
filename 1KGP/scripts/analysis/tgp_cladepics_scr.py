#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pickle

import tszip
import stdpopsim
import vcf
from liftover import get_lifter

converter = get_lifter("hg38", "hg19")
converter_back = get_lifter("hg19", "hg38")

import numpy as np
from scipy.stats import pearsonr
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

sys.path.append("..")
import dolores.simulations
import dolores.viz

col_green = "#228833"
col_red = "#EE6677"
col_purp = "#AA3377"
col_blue = "#66CCEE"
col_yellow = "#CCBB44"
col_indigo = "#4477AA"
col_grey = "#BBBBBB"


chromosome = sys.argv[1]
band = sys.argv[2]
tree_positions = sys.argv[3]
clades_to_time = sys.argv[4]
clades_to_include = sys.argv[5]
leftmargin = rightmargin = 0.12
if len(sys.argv) > 6:
    leftmargin = float(sys.argv[6])
    rightmargin = float(sys.argv[7])

name = chromosome + band + "_details"
if tree_positions != "None":
    tree_positions = [int(i) for i in tree_positions.split(sep=",")]
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

plot_off = False
full_ts = False

def get_bitset(n_set, samp_to_ind):
    e = [0] * int(len(samp_to_ind) / 2)
    for n in n_set:
        e[samp_to_ind[n]] += 1
    return e

def hg38_to_hg19(ch, pos):
    p = converter[str(ch)][int(pos)]
    if len(p) == 1:
        p = p[0][1]
    else:
        print("Cannot liftover uniquely")
        p = None
    return p

def get_chunk(ch, pos):
    chunk = None
    with open(
        "trees/relate_1000GP/trees/chunks/treeinfo_"
        + ch
        + ".txt",
        "r",
    ) as file:
        for l, line in enumerate(file):
            if l > 0:
                line = line.strip().split(";")
                if int(line[3]) <= pos < int(line[4]) and not (
                    int(line[3]) == 0 and int(line[4]) - int(line[3]) > 10000000
                ):
                    chunk = int(line[1])
    return chunk

# For paper:
# python -m tgp_cladepics_scr chr17 q21.31 43500000,43750000,43880000,44150000,44460000 995414 995414,1315245,669706,924648
# python -m tgp_cladepics_scr chr10 q22.3 81200000,81350000,81450000,81710000,82100000 971004 1000419,971004,972176,651513
# python -m tgp_cladepics_scr chr11 q11 55000000,55200000,55340000,55500000,55680000 -1 998214,684826,999618,968799 0.12 0.2
# python -m tgp_cladepics_scr chr6 p11.2 57200000,57357500,57475000,57570000,57750000 -1 1313075,686103,1004145,999726,963805 0.4 0.2
# python -m tgp_cladepics_scr chr11 q12.1 55780000,55900000,56120000,56240000,56320000 -1 1319404,691145,1000826,958808
# python -m tgp_cladepics_scr chr7 q11.21 65840000,66050000,66250000,66400000,66520000 -1 1002794,689969,1317722,973509,1001786 0.40 0.1

output_dir = "results_1000GP"
mask_loc = "files/20140520.pilot_mask.autosomes.bed"
results_file = (
    "files/clades_pvalues_all_events_annotated_lowth_varNe.csv"
)
superdups_file = "files/genomicSuperDups.txt"
svs_file = "files/common_1000g.bed"
annot_file = "files/ncbiRefSeqCurated.txt"
svgenotypes_file = "files/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"
trees_loc = "trees/relate_1000GP/trees/chunks/"
time_file = "results_1000GP/" + chromosome + band + "_time.csv"
snp_file = "results_1000GP/" + chromosome + band + "_snps.csv"
mutmap_file = "results_1000GP/" + chromosome + band + "_mutmap.csv"
cnv_file = "results_1000GP/cnvs.csv"
gwas_file = "files/gwas_catalog_v1.0.2-associations_e111_r2024-04-22.tsv"
carriers_file = "results_1000GP/" + chromosome + band + "_carriers_strict.csv"

mask = dolores.simulations.read_genome_mask_bed(mask_loc, chromosome)
(
    chromosome_lengths,
    chromosome_starts,
    centromeres,
    centromeres_mask,
) = dolores.simulations.read_genome_info(output_dir + "/chromosomes_info.txt")
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
    "ALL": 5008
}

species = stdpopsim.get_species("HomSap")
contig = species.get_contig(chromosome=chromosome, genetic_map="HapMapII_GRCh37")
recombination_map = contig.recombination_map

md, pops, groups = dolores.simulations.read_in_pop_metadata(
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


if chromosome == "chr17" or chromosome == "chr6":
    fig, axs = plt.subplots(
        12,
        1,
        figsize=(11, 12),
        gridspec_kw={"height_ratios": [0.2, 9, 2.0, 1.0, 0.5, 0.5, 1.0, 0.5, 2.0, 3.0, 1.0, 0.6]},
    )
elif chromosome == "chr10":
    fig, axs = plt.subplots(
        12,
        1,
        figsize=(11, 10),
        gridspec_kw={"height_ratios": [0.2, 4, 2.0, 1.0, 0.5, 0.5, 1.0, 0.5, 2.0, 4.0, 1.0, 0.6]},
    )
else:
    fig, axs = plt.subplots(
        12,
        1,
        figsize=(11, 12),
        gridspec_kw={
            "height_ratios": [0.2, 6, 2.0, 1.0, 0.5, 0.5, 1.0, 0.5, 2.0, 4.0, 1.0, 0.6]
        },
    )
(
    ax_scalebar,
    ax,
    ax_gene,
    ax_muts,
    ax_trees,
    ax_mask,
    ax_dups,
    ax_cnv,
    ax_1kgp,
    ax_time,
    ax_recrate,
    ax_pos,
) = axs
# plt.rcParams.update({'axes.titlesize': 'small'})

xmin = math.inf
xmax = -1
chunkmin = math.inf

if not plot_off:
    yy = 1
    ymax = 0
    bitset_all = [0] * int(len(samp_to_ind_all) / 2)
    bitset_all_strict = [0] * int(len(samp_to_ind_all) / 2)
    if os.path.exists(time_file):
        with open(results_file, "r") as file:
            for line in file:
                line = line.strip().split(sep=",")
                if line[0] == chromosome:
                    if line[21] == chromosome[3:] + band:
                        start = int(line[10])
                        end = int(line[11])
                        xmin = min(xmin, start)
                        xmax = max(xmax, end)
                        chunkmin = min(chunkmin, int(line[18]))
        print(xmin, xmax, chunkmin)
        with open(results_file, "r") as file:
            print("Reading in time file...")
            for line in file:
                line = line.strip().split(sep=",")
                if line[0] == chromosome:
                    if line[21] == chromosome[3:] + band:
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
                        print(node_id)
                        cladesize = int(line[8])
                        p1 = float(line[6])
                        p2 = float(line[7])
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
                        ts_full = tszip.decompress(
                            trees_loc
                            + chromosome
                            + "/1000GP_Phase3_mask_prene_"
                            + chromosome
                            + "_chunk"
                            + str(chunkindex)
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
                        t_ = ts_full.at_index(treeindex)
                        samples = {s for s in t.samples(node_id)}
                        bitset_true = get_bitset(samples, samp_to_ind)
                        if node_id in clades_to_include:
                            samples_ = [np.where(mapp == s)[0][0] for s in samples]
                            mrca = t_.mrca(*samples_)
                            samples__ = t_.samples(mrca)
                            for ss in samples__:
                                bitset_all[samp_to_ind_all[ss]] += 1
                            for s in samples:
                                ss = np.where(mapp == s)[0][0]
                                bitset_all_strict[samp_to_ind_all[ss]] += 1
                                # print(pop, md[ss]['ID'], ss)

                        ax.scatter([start, end], [yy+0.15, yy+0.15], marker=7, color=colors[pop], s=50)
                        ax.hlines(xmin=xmin, xmax=xmax, y=yy, color=col_grey, lw=1, zorder=-20)
                        # ax.hlines(xmin=start, xmax=left_mut, y=yy, color=col_grey, ls="dotted", lw=2)
                        # ax.hlines(xmin=right_mut, xmax=end, y=yy, color=col_grey, ls="dotted", lw=2)

                        if node_id in clades_to_include:
                            ax.text(
                                xmax+2000,
                                yy,
                                pop + ":" + str(round(cladesize / num_samples[pop], 2)) + "  ⬤",
                                ha="left",
                                va="center",
                                fontsize=7,
                            )
                        else:
                            ax.text(
                                xmax+2000,
                                yy,
                                pop + ":" + str(round(cladesize / num_samples[pop], 2)),
                                ha="left",
                                va="center",
                                fontsize=7,
                            )
                        ax.text(
                            xmin-2000,
                            yy,
                            "p1:" + str(round(p1, 1)) + ", p2:" + str(round(p2, 1)),
                            ha="right",
                            va="center",
                            fontsize=7,
                        )
                        
                        if (start < centromeres[chromosome]) and (end > centromeres[chromosome]):
                            ax.scatter(centromeres[chromosome], yy, s=50, color="black")
                            ax.hlines(
                                xmin=centromeres_mask[chromosome][0],
                                xmax=centromeres_mask[chromosome][1],
                                y=yy,
                                lw=10,
                                color="grey",
                                alpha=0.5,
                            )

                        with open(time_file, "r") as tfile:
                            for tline in tfile:
                                tline = tline.strip().split(sep=",")
                                node_id_ = int(tline[1])
                                if node_id == node_id_:
                                    mutations = [m for m in tline[23].strip().split(sep=" ")]
                                    if mutations[0] != '':
                                        mutations = set(int(m) for m in mutations)
                                    r2 = float(tline[5])
                                    left = float(tline[19])
                                    right = float(tline[20])
                                    t1 = float(tline[21])
                                    t2 = float(tline[22])
                                    if r2 > 0.90:
                                        if r2 > 0.95:
                                            ax.hlines(
                                                xmin=left,
                                                xmax=right,
                                                y=yy,
                                                color=colors[pop],
                                                lw=3,
                                            )
                                            ax.scatter(
                                                np.sort(list(mutations)),
                                                [yy] * len(mutations),
                                                marker="|",
                                                color=colors[pop],
                                                s=40,
                                            )
                                        else:
                                            ax.hlines(
                                                xmin=left,
                                                xmax=right,
                                                y=yy,
                                                color="grey",
                                                lw=3,
                                            )
                                            ax.scatter(
                                                np.sort(list(mutations)),
                                                [yy] * len(mutations),
                                                marker="|",
                                                color="grey",
                                                s=40,
                                            )
                                        if node_id in clades_to_time:
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
                        yy += 1
            pop = "ALL"
            print_all = False
            time_file_ = time_file
            if os.path.exists("results_1000GP/" + chromosome + band + "_time_ALL.csv"):
                time_file_ = "results_1000GP/" + chromosome + band + "_time_ALL.csv"
            with open(time_file_, "r") as tfile:
                for tline in tfile:
                    tline = tline.split(sep=",")
                    if tline[0] == "ALL":
                        print_all = True
                        r2 = float(tline[5])
                        left = float(tline[19])
                        right = float(tline[20])
                        t1 = float(tline[21])
                        t2 = float(tline[22])
                        cladesize = int(tline[6])
                        if right > xmax:
                            break
                        if r2 > 0.90:
                            if r2 >= 0.95:
                                c = colors[pop]
                            else:
                                c = "grey"
                            if time_all:
                                p = patches.Rectangle(
                                    (left, t1),
                                    right - left,
                                    t2 - t1,
                                    facecolor=c,
                                )
                                ymax = max(ymax, t2)
                                ax_time.add_patch(p)
                            mutations = [m for m in tline[23].strip().split(sep=" ")]
                            if mutations[0] != "":
                                mutations = set(int(m) for m in mutations)
                            if r2 > 0.95:
                                # print({m for m in mutations})
                                ax.hlines(
                                    xmin=left,
                                    xmax=right,
                                    y=yy,
                                    color=colors[pop],
                                    lw=3,
                                )
                                ax.scatter(
                                    np.sort(list(mutations)),
                                    [yy] * len(mutations),
                                    marker="|",
                                    color=colors[pop],
                                    s=40,
                                )
                            else:
                                ax.hlines(
                                    xmin=left,
                                    xmax=right,
                                    y=yy,
                                    color="grey",
                                    lw=3,
                                    ls="dotted",
                                )
                                ax.scatter(
                                    np.sort(list(mutations)),
                                    [yy] * len(mutations),
                                    marker="|",
                                    color="grey",
                                    s=40,
                                )
                if print_all:
                    ax.text(
                        xmax+2000,
                        yy,
                        pop + ":" + str(round(cladesize / num_samples[pop], 2)),
                        ha="left",
                        va="center",
                        fontsize=7,
                    )
                    ax.hlines(xmin=xmin, xmax=xmax, y=yy, color=col_grey, lw=1, zorder=-20)
                    yy += 1
            # if not print_all:
            #     ts = tszip.decompress(
            #         trees_loc
            #         + chromosome
            #         + "/1000GP_Phase3_mask_prene_"
            #         + chromosome
            #         + "_chunk"
            #         + str(chunkmin)
            #         + ".trees.tsz"
            #     )
            #     print("Getting ts at chunk", chunkmin, "at position", xmin)
            #     t = ts.at(xmin)
            #     minmaxr = math.inf
            #     tagging_SNPs = set()
            #     with open("/Users/anastasia.ignatieva/Data/results_1000GP/" + chromosome + band + "_time_ALL.csv", "a") as tfile:
            #         while t.interval[0] < xmax:
            #             maxr = 0
            #             maxr_n = -1
            #             for n in t.nodes():
            #                 if n != t.root and not t.is_sample(n):
            #                     samples_ = {s for s in t.samples(n)}
            #                     b = get_bitset(samples_, samp_to_ind_all)
            #                     if len(b) != len(bitset_all):
            #                         print(len(b), np.sum(b))
            #                         sys.exit("len(b) != len(bitset_all)")
            #                     r = pearsonr(b, bitset_all)[0]
            #                     if r > maxr:
            #                         maxr = r
            #                         maxr_n = n
            #                     if b == bitset_all:
            #                         print("Clade found at", t.interval)
            #                         break
            #             minmaxr = min(minmaxr, maxr)
            #             mutcount_clade = mutprop_clade = 0
            #             focal_mut = -1
            #             mutations = set()
            #             for m in t.mutations():
            #                 if t.is_descendant(m.node, maxr_n):
            #                     if m.node != maxr_n:
            #                         mutcount_clade += 1
            #                     else:
            #                         mutations.add(int(ts.site(m.site).position))
            #                         if focal_mut == -1:
            #                             focal_mut = int(ts.site(m.site).position)
            #             if t.num_mutations != 0:
            #                 mutprop_clade = mutcount_clade / t.num_mutations
            #             tagging_SNPs.update(mutations)
            #             counts = {g: 0 for g in groups}
            #             for s in t.samples(maxr_n):
            #                 counts[md[s]["group"]] += 1
            #             counts = {v: k / num_samples[v] for v, k in counts.items()}
            #             tbl_mut = sum([t.branch_length(s) for s in t.samples(maxr_n) if s != maxr_n])
            #             tbl = t.total_branch_length
            #             # POP,CLADE_ID,CHUNK,TREEINDEX,NODE,R2,CLADESIZE,AF_EUR,AF_AFR,AF_SAS,AF_EAS,AF_AMR,FOCAL_MUT,
            #             # MUTCOUNT_CLADE,MUTCOUNT_ALL,MUTPROP_CLADE,TBL_CLADE,TBL_ALL,TBLPROP_CLADE,
            #             # START,END,TIME_LOW,TIME_HIGH,MUTATIONS
            #             tfile.write(
            #                 "ALL"
            #                 + ","
            #                 + str(-1)
            #                 + ","
            #                 + str(chunkmin)
            #                 + ","
            #                 + str(t.index)
            #                 + ","
            #                 + str(maxr_n)
            #                 + ","
            #                 + str(maxr)
            #                 + ","
            #                 + str(t.num_samples(maxr_n))
            #                 + ","
            #                 + str(counts["EUR"])
            #                 + ","
            #                 + str(counts["AFR"])
            #                 + ","
            #                 + str(counts["SAS"])
            #                 + ","
            #                 + str(counts["EAS"])
            #                 + ","
            #                 + str(counts["AMR"])
            #                 + ","
            #                 + str(focal_mut)
            #                 + ","
            #                 + str(mutcount_clade)
            #                 + ","
            #                 + str(t.num_mutations)
            #                 + ","
            #                 + str(mutprop_clade)
            #                 + ","
            #                 + str(tbl_mut)
            #                 + ","
            #                 + str(tbl)
            #                 + ","
            #                 + str(tbl_mut / tbl)
            #                 + ","
            #                 + str(t.interval[0])
            #                 + ","
            #                 + str(t.interval[1])
            #                 + ","
            #                 + str(t.time(maxr_n))
            #                 + ","
            #                 + str(t.time(t.parent(maxr_n)))
            #                 + ","
            #                 + " ".join(str(m) for m in mutations)
            #                 + "\n"
            #             )
            #
            #             if t.index == ts.num_trees - 2:
            #                 chunkmin += 1
            #                 ts = tszip.decompress(
            #                     trees_loc
            #                     + chromosome
            #                     + "/1000GP_Phase3_mask_prene_"
            #                     + chromosome
            #                     + "_chunk"
            #                     + str(chunkmin)
            #                     + ".trees.tsz"
            #                 )
            #                 t = ts.first()
            #             t.next()
            #     print("Tagging SNPs:")
            #     with open(snp_file, "w") as sfile:
            #         for s in tagging_SNPs:
            #             sfile.write(str(s) + ",")
            #         sfile.write("\n")
            #     print("minmaxr", minmaxr)

    else:
        bitset_true_all = {}
        count = 0
        with open(results_file, "r") as file, open(time_file, "w") as tfile:
            for line in file:
                line = line.split(sep=",")
                if line[0] == chromosome:
                    if line[21][:-1] == chromosome[3:] + band:
                        pop = line[1]
                        start = int(line[10])
                        end = int(line[11])
                        chunkindex = int(line[18])
                        treeindex = int(line[19])
                        node_id = int(line[20])
                        xmin = min(xmin, start)
                        xmax = max(xmax, end)
                        chunkmin = min(chunkmin, chunkindex)
                        samp_to_ind = pickle.load(
                            open(
                                trees_loc
                                + chromosome
                                + "/1000GP_Phase3_mask_prene_"
                                + chromosome
                                + "_"
                                + pop
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
                            + "_"
                            + pop
                            + ".trees.tsz"
                        )
                        mapp = pickle.load(
                            open(
                                trees_loc
                                + chromosome
                                + "/1000GP_Phase3_mask_prene_"
                                + chromosome
                                + "_"
                                + pop
                                + "_map.pickle",
                                "rb",
                            )
                        )
                        t = ts.at_index(treeindex)
                        samples = {s for s in t.samples(node_id)}
                        bitset_true = get_bitset(samples, samp_to_ind)
                        bitset_true_all[node_id] = bitset_true
                        if node_id in clades_to_include:
                            for s in samples:
                                ss = np.where(mapp == s)[0][0]
                                bitset_all[samp_to_ind_all[ss]] += 1
        with open(results_file, "r") as file, open(time_file, "w") as tfile:
            for line in file:
                line = line.split(sep=",")
                if line[0] == chromosome:
                    if line[21][:-1] == chromosome[3:] + band:
                        pop = line[1]
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
                        bitset_true = bitset_true_all[node_id]

                        print(line[21][:-1], node_id, p1, p2, yy, cladesize)

                        ax.scatter([start, end], [yy+0.15, yy+0.15], marker=7, color=colors[pop], s=40)
                        ax.hlines(xmin=xmin, xmax=xmax, y=yy, color=col_grey, lw=1, zorder=-20)
                        # ax.hlines(
                        #     xmin=start, xmax=left_mut, y=yy, color=col_grey, ls="dotted", lw=2
                        # )
                        # ax.hlines(
                        #     xmin=right_mut, xmax=end, y=yy, color=col_grey, ls="dotted", lw=2
                        # )
                        ax.text(
                            xmax+2000,
                            yy,
                            pop + ":" + str(round(cladesize / num_samples[pop], 2)),
                            ha="left",
                            va="center",
                            fontsize=7,
                        )
                        ax.text(
                            xmin-2000,
                            yy,
                            "p1:" + str(round(p1, 1)) + ", p2:" + str(round(p2, 1)),
                            ha="right",
                            va="center",
                            fontsize=7,
                        )

                        if (start < centromeres[chromosome]) and (
                            end > centromeres[chromosome]
                        ):
                            ax.scatter(centromeres[chromosome], yy, s=50, color="black")
                            ax.hlines(
                                xmin=centromeres_mask[chromosome][0],
                                xmax=centromeres_mask[chromosome][1],
                                y=yy,
                                lw=10,
                                color="grey",
                                alpha=0.5,
                            )

                        samp_to_ind = pickle.load(
                            open(
                                trees_loc
                                + chromosome
                                + "/1000GP_Phase3_mask_prene_"
                                + chromosome
                                + "_"
                                + pop
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
                            + str(chunkmin)
                            + "_"
                            + pop
                            + ".trees.tsz"
                        )

                        t = ts.at(xmin)
                        minmaxr = math.inf

                        while t.interval[0] < xmax:
                            maxr = 0
                            maxr_n = -1
                            exact = False
                            for n in t.nodes():
                                if n != t.root and not t.is_sample(n):
                                    samples_ = {s for s in t.samples(n)}
                                    b = get_bitset(samples_, samp_to_ind)
                                    r = pearsonr(b, bitset_true)[0]
                                    if r > maxr:
                                        maxr = r
                                        maxr_n = n
                            minmaxr = min(minmaxr, maxr)
                            counts = {g: 0 for g in groups}
                            counts[pop] = t.num_samples(maxr_n)/num_samples[pop]
                            mutcount_clade = mutprop_clade = 0
                            focal_mut = -1
                            mutations = set()
                            for m in t.mutations():
                                if t.is_descendant(m.node, maxr_n):
                                    if m.node != maxr_n:
                                        mutcount_clade += 1
                                    else:
                                        mutations.add(int(ts.site(m.site).position))
                                        if focal_mut == -1:
                                            focal_mut = int(ts.site(m.site).position)
                            if t.num_mutations != 0:
                                mutprop_clade = mutcount_clade/t.num_mutations
                            tbl_mut = sum([t.branch_length(s) for s in t.nodes(maxr_n) if s != maxr_n])
                            tbl = t.total_branch_length
                            # POP,CLADE_ID,CHUNK,TREEINDEX,NODE,R2,CLADESIZE,AF_EUR,AF_AFR,AF_SAS,AF_EAS,AF_AMR,FOCAL_MUT,
                            # MUTCOUNT_CLADE,MUTCOUNT_ALL,MUTPROP_CLADE,TBL_CLADE,TBL_ALL,TBLPROP_CLADE,
                            # START,END,TIME_LOW,TIME_HIGH,MUTATIONS
                            tfile.write(
                                pop
                                + ","
                                + str(node_id)
                                + ","
                                + str(chunkmin)
                                + ","
                                + str(t.index)
                                + ","
                                + str(maxr_n)
                                + ","
                                + str(maxr)
                                + ","
                                + str(t.num_samples(maxr_n))
                                + ","
                                + str(counts["EUR"])
                                + ","
                                + str(counts["AFR"])
                                + ","
                                + str(counts["SAS"])
                                + ","
                                + str(counts["EAS"])
                                + ","
                                + str(counts["AMR"])
                                + ","
                                + str(focal_mut)
                                + ","
                                + str(mutcount_clade)
                                + ","
                                + str(t.num_mutations)
                                + ","
                                + str(mutprop_clade)
                                + ","
                                + str(tbl_mut)
                                + ","
                                + str(tbl)
                                + ","
                                + str(tbl_mut/tbl)
                                +","
                                + str(t.interval[0])
                                + ","
                                + str(t.interval[1])
                                + ","
                                + str(t.time(maxr_n))
                                + ","
                                + str(t.time(t.parent(maxr_n)))
                                + ","
                                + " ".join(str(m) for m in mutations)
                                + "\n"
                            )
                            if maxr > 0.90:
                                if maxr > 0.95:
                                    ax.hlines(
                                        xmin=t.interval[0],
                                        xmax=t.interval[1],
                                        y=yy,
                                        color=colors[pop],
                                        lw=3,
                                    )
                                    ax.scatter(
                                        np.sort(list(mutations)),
                                        [yy] * len(mutations),
                                        marker="|",
                                        color=colors[pop],
                                        s=40,
                                    )
                                else:
                                    ax.hlines(
                                        xmin=t.interval[0],
                                        xmax=t.interval[1],
                                        y=yy,
                                        color="grey",
                                        lw=3,
                                    )
                                    ax.scatter(
                                        np.sort(list(mutations)),
                                        [yy] * len(mutations),
                                        marker="|",
                                        color="grey",
                                        s=40,
                                    )
                                if node_id in clades_to_time:
                                    if maxr >= 0.95:
                                        c = colors[pop]
                                    else:
                                        c = "grey"
                                    p = patches.Rectangle(
                                        (t.interval[0], t.time(maxr_n)),
                                        t.interval[1] - t.interval[0],
                                        t.time(t.parent(maxr_n)) - t.time(maxr_n),
                                        facecolor=c,
                                    )
                                    ymax = max(ymax, t.time(t.parent(maxr_n)))
                                    ax_time.add_patch(p)
                                # if exact and cladesize < 500 and len(mutations) > 0:
                                #     print(mutations)
                            if t.index == ts.num_trees - 2:
                                ts = tszip.decompress(
                                    trees_loc
                                    + chromosome
                                    + "/1000GP_Phase3_mask_prene_"
                                    + chromosome
                                    + "_chunk"
                                    + str(chunkmin + 1)
                                    + "_"
                                    + pop
                                    + ".trees.tsz"
                                )
                                t = ts.first()
                            t.next()
                        print("minmaxr", minmaxr)
                        yy += 1
        print("Total samples in clades:", np.sum(bitset_all))
        print("Total individuals:", len([b for b in bitset_all if b > 0]))

        if full_ts:
            ts = tszip.decompress(
                trees_loc
                + chromosome
                + "/1000GP_Phase3_mask_prene_"
                + chromosome
                + "_chunk"
                + str(chunkmin)
                + ".trees.tsz"
            )
            print("Getting ts at chunk", chunkmin, "at position", xmin)
            t = ts.at(xmin)
            minmaxr = math.inf
            with open(time_file, "a") as tfile:
                while t.interval[0] < xmax:
                    maxr = 0
                    maxr_n = -1
                    for n in t.nodes():
                        if n != t.root and not t.is_sample(n):
                            samples_ = {s for s in t.samples(n)}
                            b = get_bitset(samples_, samp_to_ind_all)
                            if len(b) != len(bitset_all):
                                print(len(b), np.sum(b))
                                sys.exit("len(b) != len(bitset_all)")
                            r = pearsonr(b, bitset_all)[0]
                            if r > maxr:
                                maxr = r
                                maxr_n = n
                            if b == bitset_all:
                                print("Clade found at", t.interval)
                                break
                    minmaxr = min(minmaxr, maxr)
                    mutcount_clade = mutprop_clade = 0
                    focal_mut = -1
                    mutations = set()
                    for m in t.mutations():
                        if t.is_descendant(m.node, maxr_n):
                            if m.node != maxr_n:
                                mutcount_clade += 1
                            else:
                                mutations.add(int(ts.site(m.site).position))
                                if focal_mut == -1:
                                    focal_mut = int(ts.site(m.site).position)
                    if t.num_mutations != 0:
                        mutprop_clade = mutcount_clade / t.num_mutations
                    counts = {g: 0 for g in groups}
                    for s in t.samples(maxr_n):
                        counts[md[s]["group"]] += 1
                    counts = {v: k / num_samples[v] for v, k in counts.items()}
                    tbl_mut = sum([t.branch_length(s) for s in t.nodes(maxr_n) if s != maxr_n])
                    tbl = t.total_branch_length
                    tfile.write(
                        "ALL"
                        + ","
                        + str(-1)
                        + ","
                        + str(chunkmin)
                        + ","
                        + str(t.index)
                        + ","
                        + str(maxr_n)
                        + ","
                        + str(maxr)
                        + ","
                        + str(t.num_samples(maxr_n))
                        + ","
                        + str(counts["EUR"])
                        + ","
                        + str(counts["AFR"])
                        + ","
                        + str(counts["SAS"])
                        + ","
                        + str(counts["EAS"])
                        + ","
                        + str(counts["AMR"])
                        + ","
                        + str(focal_mut)
                        +","
                        + str(mutcount_clade)
                        + ","
                        + str(t.num_mutations)
                        + ","
                        + str(mutprop_clade)
                        + ","
                        + str(tbl_mut)
                        + ","
                        + str(tbl)
                        + ","
                        + str(tbl_mut/tbl)
                        + ","
                        + str(t.interval[0])
                        + ","
                        + str(t.interval[1])
                        + ","
                        + str(t.time(maxr_n))
                        + ","
                        + str(t.time(t.parent(maxr_n)))
                        + ","
                        + " ".join(str(m) for m in mutations)
                        + "\n"
                    )

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
            print("minmaxr", minmaxr)

with open(carriers_file, "w") as cfile:
    for i in inds:
        b = bitset_all_strict[inds[i]]
        for k in md:
            if md[k]['ID'] == i:
                pop1 = md[k]['group']
                pop2 = md[k]['population']
                break
        cfile.write(i + "," + pop1 + "," + pop2 + "," + str(b) + "\n")

r = xmax - xmin
xmin = xmin - leftmargin * r
xmax = xmax + rightmargin * r
print(xmin, xmax)

xticks = [i*100000 for i in range(0, 3000) if xmin <= i*100000 <= xmax]
xlabels = [str(i/1000000) for i in xticks]
plt.setp(axs, xticks=xticks, xticklabels=xlabels)

ax.set_title("DoLoReS: Genomic span of clades, " + chromosome[3:] + band)
ax.set_xlim([xmin, xmax])
ax.tick_params(axis="x", which="both", labelbottom=False)
ax.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

"""
Age estimate of inversion
"""

ax_time.set_title("Age estimate")
ax_time.set_xlim([xmin, xmax])
ax_time.tick_params(axis="x", which="both", labelbottom=False)
if plot_off:
    ymax = 10000
ax_time.set_ylim([0, ymax*1.05])

"""
Segmental duplications >1000bp from UCSC browser
"""

if not plot_off:
    dups = dolores.simulations.read_superdups(
        superdups_file,
        xmin,
        xmax,
        chromosome,
    )
    for dup in dups:
        if dup[2] in [1, 2, 3]:
            print(dup, dup[1]-dup[0])
        if dup[2] == 1:
            p = patches.Rectangle(
                (dup[0], 0.5),
                dup[1] - dup[0],
                0.4,
                linewidth=0,
                facecolor=dup[3],
                hatch="/////",
                edgecolor="white",
            )
            ax_dups.add_patch(p)
        elif dup[2] == 2:
            # inverted repeats
            p = patches.Rectangle(
                (dup[0], 0),
                dup[1] - dup[0],
                0.4,
                linewidth=0,
                facecolor=dup[3],
                hatch="/////",
                edgecolor="white",
            )
            ax_dups.add_patch(p)
        elif dup[2] == 3:
            # inverted repeats
            p = patches.Rectangle(
                (dup[0], 0),
                dup[1] - dup[0],
                0.4,
                linewidth=0,
                facecolor=dup[3],
                hatch="\\\\\\\\\\",
                edgecolor="white",
            )
            ax_dups.add_patch(p)
ax_dups.set_ylim([-0.1, 1.0])
ax_dups.set_xlim([xmin, xmax])
ax_dups.set_title("Segmental duplications (>1kb, >90% similarity)")
ax_dups.text(
    x=xmin + 10000,
    y=0.2,
    s="INVERTED",
    ha="left",
    va="center",
    fontsize=7,
)
ax_dups.text(
    x=xmin + 10000,
    y=0.7,
    s="DIRECT",
    ha="left",
    va="center",
    fontsize=7,
)
ax_dups.tick_params(axis="x", which="both", labelbottom=False)
ax_dups.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

"""
SVs in 1000GP call set
"""

if not plot_off:
    svs = dolores.simulations.read_1KGP_svs(
        svs_file,
        xmin,
        xmax,
        chromosome,
    )
    print("SVs", svs)
    yy = 1.5
    svs_ = []
    svs__ = []
    for sv in svs:
        if sv[1] - sv[0] + 1 > 100 and (sv[0], sv[1], sv[2]) not in svs_:
            svs_.append((sv[0], sv[1], sv[2]))
            if yy == 0.0:
                yy = 1.5
            else:
                yy = yy - 0.5
            if sv[2] == "deletion":
                s = "DEL"
            elif sv[2] == "duplication":
                s = "DUP"
            elif sv[2] == "copy":
                s = "CNV"
            elif sv[2] == "alu":
                s = "ALU"
            elif sv[2] == "line1":
                s = "LINE1"
            elif sv[2] == "sva":
                s = "SVA"
            else:
                s = "OTHER"
            svs__.append((sv[0], sv[1], s, yy))
            p = patches.Rectangle(
                (sv[0], yy),
                sv[1] - sv[0],
                0.4,
                linewidth=0,
                facecolor="black",
            )
            ax_1kgp.add_patch(p)
            if sv[1] < xmax:
                ax_1kgp.text(
                    x=sv[1],
                    y=yy + 0.2,
                    s=s,
                    ha="left",
                    va="center",
                    fontsize=7,
                )
ax_1kgp.set_title("1KGP SV call set (>100bp)")
ax_1kgp.set_ylim([-0.1, 2.0])
ax_1kgp.set_xlim([xmin, xmax])
ax_1kgp.tick_params(axis="x", which="both", labelbottom=False)
ax_1kgp.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

vcf_reader_sv = vcf.Reader(
    filename=svgenotypes_file,
    compressed=True,
)
records_sv = vcf_reader_sv.fetch(int(chromosome[3:]), int(xmin), int(xmax))
for record2 in records_sv:
    bitset_sv = [0] * int(len(samp_to_ind_all) / 2)
    M = 0
    for sample in record2.samples:
        i = inds[sample.sample]
        bitset_sv[i] += int(sample["GT"][0])
        bitset_sv[i] += int(sample["GT"][2])
        M = M + int(sample["GT"][0]) + int(sample["GT"][2])
    r = pearsonr(bitset_sv, bitset_all)[0]
    # inx = [i for i in range(len(bitset_all)) if bitset_all[i] != 0]
    # bitset_sv_ = [bitset_sv[i] for i in inx]
    # bitset_all_ = [bitset_all[i] for i in inx]
    # r_ = pearsonr(bitset_sv_, bitset_all_)[0]
    for p1, p2, t, yy in svs__:
        if abs(record2.POS - p1) <= 5 and t == record2.var_subtype:
            ax_1kgp.text(
                x=p2,
                y=yy + 0.2,
                s=t + ":" + "%.2f" % round(r**2, 2),
                ha="left",
                va="center",
                fontsize=7,
            )
    print(record2.POS, record2.var_subtype, M, r**2)

"""
GWAS hits
"""
# nlines = 0
# with open(gwas_file, "r") as file:
#     for line in file:
#         nlines += 1
# with tqdm.tqdm(total=nlines) as pbar:
#     with open(gwas_file, "r") as file:
#         for l, line in enumerate(file):
#             if l > 0:
#                 line = line.strip().split(sep="\t")
#                 if line[11] != "" and line[11] != "X" and line[11] != "Y":
#                     if ("x" not in line[11]) and (";" not in line[11]):
#                         ch = int(line[11])
#                         if "chr" + str(ch) == chromosome:
#                             trait = line[7]
#                             pos = int(line[12])
#                             pos_ = hg38_to_hg19(ch, pos)
#                             if pos_ is not None:
#                                 if xmin - 1000000 <= pos_ <= xmax + 1000000:
#                                     # print(trait, ch, pos, pos_)
#                                     chunk = get_chunk(chromosome, pos_)
#                                     ts = tszip.decompress(
#                                         trees_loc
#                                         + chromosome
#                                         + "/1000GP_Phase3_mask_prene_"
#                                         + chromosome
#                                         + "_chunk"
#                                         + str(chunk)
#                                         + ".trees.tsz"
#                                     )
#                                     t = ts.at(pos_)
#                                     for m in t.mutations():
#                                         if ts.site(m.site).position == pos_:
#                                             samples = {s for s in t.samples(m.node)}
#                                             bitset = get_bitset(samples, samp_to_ind_all)
#                                             r = pearsonr(bitset, bitset_all)[0]
#                                             r2 = r**2
#                                             if r2 > 0.6:
#                                                 print("===>", r2, trait, ch, pos, pos_, line[23], "<===")
#             pbar.update(1)

"""
RefSeq genes
"""
if not plot_off:
    genes = set()
    yy = 1.5
    with open(annot_file, "r") as file:
        for line in file:
            line = line.split()
            ch = line[2]
            if ch == chromosome:
                gene = line[12]
                start = int(line[4])
                end = int(line[5])
                if start < xmax and end > xmin:
                    if gene not in genes:
                        if yy == 0.0:
                            yy = 1.5
                        else:
                            yy = yy - 0.5
                        genes.add(gene)
                        if end < xmax - 50000:
                            ax_gene.text(
                                x=end,
                                y=yy + 0.15,
                                s=gene,
                                ha="left",
                                va="center",
                                fontsize=7,
                            )
                    rect = patches.Rectangle(
                        (start, yy),
                        end - start,
                        0.3,
                        linewidth=None,
                        edgecolor=None,
                        facecolor="black",
                    )
                    ax_gene.add_patch(rect)

ax_gene.set_title("NCBI RefSeq genes")
ax_gene.set_ylim([-0.1, 2.1])
ax_gene.set_xlim([xmin, xmax])
ax_gene.tick_params(axis="x", which="both", labelbottom=False)
ax_gene.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

"""
Genomic mask
"""
if not plot_off:
    masked = []
    for start, end in mask:
        if start < xmax and end > xmin:
            rect = patches.Rectangle(
                (start, 0),
                end - start,
                1,
                linewidth=None,
                edgecolor=None,
                facecolor="black",
                zorder=-10,
            )
            ax_mask.add_patch(rect)
            masked = masked + [p for p in range(start, end + 1)]
ax_mask.set_ylim([0, 1])
ax_mask.set_xlim([xmin, xmax])
ax_mask.set_title("Masked regions")
ax_mask.tick_params(axis="x", which="both", labelbottom=False)
ax_mask.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

"""
Mutations not mapping
"""
if not plot_off:
    if os.path.exists(mutmap_file):
        moving_averages = []
        moving_averages_x = []
        with open(mutmap_file, "r") as file:
            for line in file:
                line = line.split(",")
                moving_averages_x.append(float(line[0]))
                moving_averages.append(float(line[1]))
    else:
        with open(mutmap_file, "w") as file:
            filename = "trees/relate_1000GP/ancmut/1000GP_Phase3_mask_prene_" + chromosome
            if os.path.exists(filename + ".mut.gz"):
                muts_all, muts = dolores.simulations.read_relate_muts(filename + ".mut.gz", xmin, xmax, masked)
            else:
                muts = []
                muts_all = []
                for p in range(3):
                    if os.path.exists(filename + "_part" + str(p) + ".mut.gz"):
                        M_all, M = dolores.simulations.read_relate_muts(filename + "_part" + str(p) + ".mut.gz", xmin, xmax, masked)
                        muts = muts + M
                        muts_all = muts_all + M_all
            window_size = 10000
            i = 0
            moving_averages = []
            moving_averages_x = []
            muts_ = [0]*(int(xmax - xmin) + 1)
            muts_all_ = [0]*(int(xmax - xmin) + 1)
            for m in muts:
                muts_[m-int(xmin)] = 1
            for m in muts_all:
                muts_all_[m-int(xmin)] = 1
            while i < len(muts_) - window_size + 1:
                window = muts_[i : i + window_size]
                window_ = muts_all_[i : i + window_size]
                if sum(window_) < 10:
                    window_average = 0
                else:
                    window_average = round(sum(window) / sum(window_), 2)
                moving_averages.append(window_average)
                moving_averages_x.append(xmin + i + window_size/2)
                file.write(
                    str(xmin + i + window_size/2)
                    + ","
                    + str(window_average)
                    + "\n"
                )
                i += 1
    ax_muts.plot(moving_averages_x, moving_averages, color="black")
ax_muts.set_title("Proportion of mutations not mapping to local tree")
ax_muts.set_ylim([0, 1])
ax_muts.set_xlim([xmin, xmax])
ax_muts.tick_params(axis="x", which="both", labelbottom=False)

"""
Recombination rate map
"""

ratemap_ = recombination_map.slice(left=xmin, right=xmax)
step = 2000
ratex = np.linspace(xmin, xmax, int((xmax - xmin)/step))
ratey = [(recombination_map.get_cumulative_mass(ratex[i]) - recombination_map.get_cumulative_mass(ratex[i-1]))/(ratex[i] - ratex[i-1]) for i in range(1, len(ratex))]
# ax_recrate.stairs(ratemap_.rate, ratemap_.position, color="black")
ax_recrate.stairs(ratey, ratex, color="black")
ax_recrate.set_xlim([xmin, xmax])
ax_recrate.set_xlabel("Genome position (Mb) on " + chromosome)
ax_recrate.set_title("Recombination rate")

"""
Tree breaks
"""

if not plot_off:
    with open(trees_loc + "treeinfo_" + chromosome + ".txt", "r") as file:
        for i, line in enumerate(file):
            if i > 0:
                line = line.split(sep=";")
                p = int(line[4])
                if xmin <= p <= xmax:
                    ax_trees.vlines(ymin=0, ymax=1, x=p, lw=0.5, color="black")
ax_trees.set_ylim([0, 1])
ax_trees.set_xlim([xmin, xmax])
ax_trees.set_title("Local tree breakpoints")
ax_trees.tick_params(axis="x", which="both", labelbottom=False)
ax_trees.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

"""
CNVs
"""

# positions = {
#     "EUR": 0.8,
#     "AFR": 0.6,
#     "EAS": 0.4,
#     "SAS": 0.2,
#     "AMR": 0.0,
# }
if not plot_off:
    with open(cnv_file, "r") as file:
        for line in file:
            line = line.split(sep=",")
            if line[0] == chromosome:
                p1 = int(line[3])
                p2 = int(line[4])
                pop = line[1]
                if p1 < xmax and p2 > xmin:
                    # rect = patches.Rectangle(
                    #     (p1, positions[pop]),
                    #     p2-p1,
                    #     0.15,
                    #     linewidth=None,
                    #     edgecolor=None,
                    #     facecolor=colors[pop],
                    # )
                    # ax_cnv.add_patch(rect)
                    ax_cnv.scatter((p1+p2)/2, 0, marker="s", s=8, color="black")
ax_cnv.set_ylim([-0.5,0.5])
ax_cnv.set_xlim([xmin, xmax])
ax_cnv.set_title("DoLoReS: Inferred CNVs")
ax_cnv.tick_params(axis="x", which="both", labelbottom=False)
ax_cnv.tick_params(
    axis="y",
    which="both",
    labelleft=False,
    left=False,
)

"""
50kb length marker
"""

# ax_scalebar.set_title(
#     chromosome + " (" + str(int(xmin)) + ", " + str(int(xmax)) + ")",
# )
rect = patches.Rectangle(
    (xmin + (xmax - xmin) * 0.05, 0.0),
    50000,
    1.0,
    linewidth=None,
    edgecolor=None,
    facecolor="black",
)
ax_scalebar.add_patch(rect)
ax_scalebar.text(
    x=xmin + (xmax - xmin) * 0.05 + 52000, y=0.5, s="50kb", ha="left", va="center"
)
ax_scalebar.set_ylim([0, 1])
ax_scalebar.set_xlim([xmin, xmax])
ax_scalebar.axis("off")

"""
Points at which trees are plotted
"""

if tree_positions != "None":
    ax_pos.scatter(
        tree_positions, [0] * len(tree_positions), s=10, color="black", marker="s"
    )
    ax_pos.set_ylim([-0.2, 2.0])
    ax_pos.set_xlim([xmin, xmax])
ax_pos.axis("off")

alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]
for i, axx in enumerate(axs):
    if i > 0:
        axx.annotate(alphabet[i-1], xy=(-0.06,1.06), xytext=(0,0), xycoords=('axes fraction', 'axes fraction'),
                    textcoords='offset points', size=14, ha="left", va="top", weight="bold")

plt.tight_layout()
plt.subplots_adjust(hspace=0.5)
plt.savefig(output_dir + "/" + name + ".png", dpi=300, bbox_inches="tight", pad_inches=0.1)
