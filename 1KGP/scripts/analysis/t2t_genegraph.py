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
import matplotlib.patches as patches
from matplotlib.lines import Line2D

sys.path.append("..")

col_green = "#228833"
col_red = "#EE6677"
col_purp = "#AA3377"
col_blue = "#66CCEE"
col_yellow = "#CCBB44"
col_indigo = "#4477AA"
col_grey = "#BBBBBB"

"""
python -m t2t_genegraph
"""

ann_dir = "t2t/annotations/"
output_dir = "t2t/plots/"
# genes_to_label = {
#     "ZMIZ1",
#     "PPIF",
#     "ZCCHC24",
#     "EIF5A",
#     "EIF5AL1",
#     "SFTPA2",
#     "MBL3P",
#     "SFTPA1",
#     "LINC02679",
#     "NUTM2B",
#     "NUTM2B-AS1",
#     "NUTM2E",
#     "CTSLP6",
#     "PGGT1BP2",
#     "BMS1P21",
#     "MBL1P",
#     "SFTPD",
#     "SFTPD-AS1",
#     "TMEM254-AS1",
#     "TMEM254",
#     "PLAC9",
#     "ANXA11",
#     "LINC00857",
#     "MAT1A",
#     "ZNF519P1",
#     "DYDC1",
#     "DYDC2",
#     "PRXL2A",
#     "TSPAN14",
# }
# col_pal = [
#     "ZMIZ1",
#     "PPIF",
#     "ZCCHC24",
#     "EIF5A",
#     "EIF5AL1",
#     "SFTPA2",
#     "SFTPA1",
#     "LINC02679",
#     "NUTM2B",
#     "NUTM2B-AS1",
#     "NUTM2E",
#     "SFTPD",
#     "SFTPD-AS1",
#     "TMEM254-AS1",
#     "TMEM254",
#     "PLAC9",
#     "ANXA11",
#     "LINC00857",
#     "MAT1A",
#     "DYDC1",
#     "DYDC2",
#     "PRXL2A",
#     "TSPAN14",
# ]
# clrs = sns.color_palette("husl", n_colors=len(col_pal))
# col_pal = {col_pal[i]: clrs[i] for i in range(len(col_pal))}

band = sys.argv[1]
startgene = sys.argv[2]
endgene = sys.argv[3]
name = str(band) + "_genegraph_" + startgene + "_" + endgene + "_updated.pdf"
genes = {startgene, endgene}
startgene_dir = endgene_dir = None
genes_to_label = []

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

with open("files/counts_HPRC_updated.csv", "r") as file:
    i = 0
    IDs = {}
    carriers = {}
    for l, line in enumerate(file):
        line = line.strip().split(",")
        if l == 0:
            for ll in line:
                if ll != "":
                    IDs[ll] = i
                    i += 1
        else:
            carriers[line[0]] = [ll for ll in line[1:] if ll != ""]

fig, axs = plt.subplots(
    1,
    1,
    figsize=(12, 20),
    gridspec_kw={
        "height_ratios": [1.0]
    },
)

directory = os.fsencode(ann_dir)
filelist = glob.glob(os.path.join(ann_dir, '*.gtf.gz'))
filecount = 0
infofile = output_dir + band + "_genegraph_" + startgene + "_" + endgene + ".csv"
if not os.path.exists(infofile):
    with open(infofile, "w") as ofile:
        # for file in os.listdir(directory):
        for f, filename in enumerate(sorted(filelist)):
            search = False
            print(filename)
            coords = {}
            with gzip.open(filename, "rt") as afile:
                for l, line in enumerate(afile):
                    line = line.strip()
                    if l == 0:
                        line = line.split(" ")
                        ref = line[1]
                    if l > 4:
                        line = line.split("\t")
                        if line[2] == "gene":
                            r = line[0]
                            s = int(line[3])
                            e = int(line[4])
                            st = line[6]
                            details_ = line[-1].split(";")
                            details = {}
                            for d in details_:
                                d = [dd.strip() for dd in d.split("\"")]
                                if len(d) > 1 and d[0] not in details:
                                    details[d[0]] = d[1]
                            if "gene_name" in details:
                                gene_name = details["gene_name"]
                                if gene_name in genes:
                                    print(gene_name, r)
                                    if gene_name in coords:
                                        print("Gene on two contigs", coords[gene_name][0], r)
                                    else:
                                        coords[gene_name] = [r, s, e, st]
                                    if f == 0:
                                        if gene_name == startgene:
                                            startgene_dir = st
                                        else:
                                            endgene_dir = st

            if coords[startgene][0] == coords[endgene][0]:
                # Same contig
                read = coords[startgene][0]
                start = min(coords[startgene][1], coords[endgene][1])
                end = max(coords[startgene][2], coords[endgene][2])
                strand = coords[startgene][3]
                rev = coords[startgene][1] > coords[endgene][1]
                # print(read, start, end)
                with gzip.open(filename, "rt") as afile:
                    for l, line in enumerate(afile):
                        if l > 4:
                            line = line.strip()
                            line = line.split("\t")
                            if line[2] == "gene":
                                r = line[0]
                                if r == read:
                                    s = int(line[3])
                                    e = int(line[4])
                                    st = line[6]
                                    if e >= start and s <= end:
                                        details_ = line[-1].split(";")
                                        details = {}
                                        for d in details_:
                                            d = [dd.strip() for dd in d.split('"')]
                                            if len(d) > 1 and d[0] not in details:
                                                details[d[0]] = d[1]
                                        if "gene_name" in details:
                                            gene_name = details["gene_name"]
                                            if gene_name in coords and gene_name not in genes:
                                                print("Already in coords")
                                                print(gene_name)
                                                print(coords[gene_name])
                                                print([read, s, e, st])
                                            else:
                                                coords[gene_name] = [read, s, e, st]
                # print(coords)
                if rev:
                    start = -start
                    end = -end
                for gene_name, (r, s, e, st) in coords.items():
                    ss = s
                    ee = e
                    if rev:
                        s = -s
                        e = -e
                        s = s - end
                        e = e - end
                    else:
                        s = s - start
                        e = e - start

                    ofile.write(
                        ref
                        + ","
                        + gene_name
                        + ","
                        + str(r)
                        + ","
                        + str(ss)
                        + ","
                        + str(s)
                        + ","
                        + str(ee)
                        + ","
                        + str(e)
                        + ","
                        + str(st)
                        + ","
                        + str(rev)
                        + ","
                        + "same"
                        + "\n"
                    )
            else:
                print("----> Different contigs <----")
                read1 = coords[startgene][0]
                read2 = coords[endgene][0]
                for i, read, gene, gene_dir in [(0, read1, startgene, startgene_dir), (1, read2, endgene, endgene_dir)]:
                    print("Checking contig", read)
                    if (i == 0 and coords[gene][3] == gene_dir) or (i == 1 and coords[gene][3] != gene_dir):
                        start = coords[gene][1]
                        end = 1000000000
                        if i == 0:
                            add = None
                            rev = False
                        else:
                            rev = True
                    else:
                        start = 0
                        end = coords[gene][2]
                        if i == 0:
                            rev = True
                            add = coords[gene][2]
                        else:
                            rev = False
                    strand = coords[gene][3]
                    end_ = 0
                    print(read, start, end, strand)
                    with gzip.open(filename, "rt") as afile:
                        for l, line in enumerate(afile):
                            if l > 4:
                                line = line.strip()
                                line = line.split("\t")
                                if line[2] == "gene":
                                    r = line[0]
                                    if r == read:
                                        s = int(line[3])
                                        e = int(line[4])
                                        st = line[6]
                                        if e >= start and s <= end:
                                            details_ = line[-1].split(";")
                                            details = {}
                                            for d in details_:
                                                d = [dd.strip() for dd in d.split('"')]
                                                if len(d) > 1 and d[0] not in details:
                                                    details[d[0]] = d[1]
                                            if "gene_name" in details:
                                                gene_name = details["gene_name"]
                                                if gene_name in coords and gene_name not in genes:
                                                    print("Already in coords")
                                                    print(coords[gene_name])
                                                    print([read, s, e, st])
                                                else:
                                                    end_ = max(end_, e)
                                                    coords[gene_name] = [read, s, e, st]
                    # print(coords)
                    end = end_
                    if add is None:
                        add = end_ - start
                    if rev:
                        start = -start
                        end = -end
                    for gene_name, (r, s, e, st) in coords.items():
                        if r == read:
                            ss = s
                            ee = e
                            if rev:
                                s = -s
                                e = -e
                                s = s - end
                                e = e - end
                            else:
                                s = s - start
                                e = e - start

                            ofile.write(
                                ref
                                + ","
                                + gene_name
                                + ","
                                + str(r)
                                + ","
                                + str(ss)
                                + ","
                                + str(s + i * add)
                                + ","
                                + str(ee)
                                + ","
                                + str(e + i * add)
                                + ","
                                + str(st)
                                + ","
                                + str(rev)
                                + ","
                                + "split"
                                + "\n"
                            )
            print(filecount)
            filecount += 1

axs.hlines(
    xmin=p1,
    xmax=p2,
    y=1,
    linewidth=5,
    color="black",
)

yy = 1
prevref = None
inds = []
with open(infofile, "r") as ofile:
    for l, line in enumerate(ofile):
        line = line.strip().split(",")
        if l == 0:
            ref_ = line[0]
        ref = line[0]
        s = int(line[4])
        if ref == ref_:
            gene_name = line[1]
            genes_to_label.append(gene_name)
            inds.append(s)
genes_to_label = [x for _, x in sorted(zip(inds, genes_to_label))]
print(genes_to_label)
clrs = sns.color_palette("husl", n_colors=len(genes_to_label))
col_pal = {genes_to_label[i]: clrs[i] for i in range(len(genes_to_label))}
with open(infofile, "r") as ofile:
    for line in ofile:
        line = line.strip().split(",")
        ref = line[0]
        gene_name = line[1]
        r = line[2]
        s = int(line[4])
        e = int(line[6])
        st = line[7]
        rev = line[8]
        if rev == "False" or rev == "FALSE":
            rev = False
        else:
            rev = True
        if prevref != ref:
            yy -= 1
            prevref = ref
            reflabel = ref.split(".")
            if len(reflabel) >= 3:
                n = reflabel[0] + "." + reflabel[2]
                if n in IDs:
                    reflabel_ = n + " " + str(carriers[band][IDs[n]]) + " "
                else:
                    reflabel_ = n + " -- "
            else:
                reflabel_ = reflabel[0] + " "
            axs.text(
                1,
                yy,
                reflabel_,
                ha="right",
                va="center",
                fontsize=10,
            )

        if (st == "+" and rev) or (st == "-" and not rev):
            sym = ">"
        else:
            sym = "<"

        col = "white"
        if gene_name in col_pal:
            col = col_pal[gene_name]
        axs.hlines(
            xmin=s,
            xmax=e,
            y=yy,
            color="black",
            lw=1,
        )
        axs.hlines(
            xmin=s,
            xmax=e,
            y=yy,
            color=col,
            lw=8,
            zorder=-20,
        )
        axs.text(
            (s + e) / 2,
            yy,
            sym,
            ha = "center",
            va = "center",
            fontsize = 12,
            )
        for x in [s, e]:
            axs.text(
                x,
                yy,
                "|",
                ha="center",
                va="center",
                fontsize=12,
            )
        if yy == 0 and gene_name in genes_to_label and abs(e-s) > 1000:
            axs.text(
                (s+e)/2,
                1.5,
                gene_name,
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
                )

plt.axis("off")
plt.tight_layout()
plt.savefig(output_dir + "/" + name, dpi=300, bbox_inches="tight")




