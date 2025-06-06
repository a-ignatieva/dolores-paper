#!/usr/bin/env python
# coding: utf-8

import math

output_dir = "results_1000GP"
results_file = (
    "files/clades_pvalues_all_events_annotated_varNe.csv"
)
genes_file = (
    "files/genes.csv"
)
annot_file = "files/ncbiRefSeqCurated.txt"

bands = {}
with open(results_file, "r") as file:
    for l, line in enumerate(file):
        if l > 0:
            line = line.split(sep=",")
            chromosome = line[0]
            band = line[21][:-1]
            if band not in bands:
                bands[band] = [math.inf, 0, chromosome]
            start = int(line[10])
            end = int(line[11])
            bands[band][0] = min(bands[band][0], start)
            bands[band][1] = max(bands[band][1], end)
for band in bands:
    print(band, bands[band], bands[band][1] - bands[band][0])

"""
RefSeq genes
"""

with open(annot_file, "r") as file, open(genes_file, "w") as out:
    out.write("chromosome,band,gene,clades_start,clades_end,clades_span,gene_start,gene_end,gene_span,prop_overlap_of_clades,prop_overlap_of_gene,\n")
    for line in file:
        line = line.split()
        ch = line[2]
        gene = line[12]
        start = int(line[4])
        end = int(line[5])
        for band in bands:
            if bands[band][2] == ch:
                xmin = bands[band][0]
                xmax = bands[band][1]
                a1 = xmax - max(xmin, start)
                a2 = xmax - min(xmax, end)
                A = (a1 - a2)/(xmax - xmin)
                B = (a1 - a2)/(end - start)
                if start < xmax and end > xmin:
                    out.write(
                        ch + ","
                        + band + ","
                        + gene + ","
                        + str(xmin) + ","
                        + str(xmax) + ","
                        + str(xmax - xmin) + ","
                        + str(start) + ","
                        + str(end) + ","
                        + str(end - start) + ","
                        + str(round(A, 5)) + ","
                        + str(round(B, 5)) + "\n"
                    )
