#!/usr/bin/env python
# coding: utf-8

import sys

import tszip
import math

import numpy as np

sys.path.append("..")
import argbd.simulations

band = sys.argv[1]
if "p" in band:
    chromosome = "chr" + band.split("p")[0]
else:
    chromosome = "chr" + band.split("q")[0]
print(chromosome, band)

output_dir = "results_1000GP"
results_file = "files/clades_pvalues_all_events_annotated_varNe.csv"
trees_loc = "trees/relate_1000GP/trees/chunks/"
time_file = (
    "results_1000GP/" + "chr" + band + "_time.csv"
)
freq_file = (
    "results_1000GP/" + "chr" + band + "_freq.csv"
)
info_file = (
    "results_1000GP/clades_info_all.csv"
)
superdups_file = "files/genomicSuperDups.txt"
cnv_file = "files/cnvs.csv"
(
    chromosome_lengths,
    chromosome_starts,
    centromeres,
    centromeres_mask,
) = argbd.simulations.read_genome_info(output_dir + "/chromosomes_info.txt")

num_samples = {
    "EUR": 1006,
    "AFR": 1322,
    "EAS": 1008,
    "SAS": 978,
    "AMR": 694,
    "ASIA": 1008,
    "ADMIX": 694,
    "ALL": 5008,
}
populations = ["EUR", "AFR", "SAS", "EAS", "AMR"]

md, pops, groups = argbd.simulations.read_in_pop_metadata(
    "trees/relate_1000GP/data/1000GP_Phase3.poplabels"
)

xmin = math.inf
xmax = -1

with open(results_file, "r") as file, open(freq_file, "w") as ffile:
    print("Reading in time file...")
    for line in file:
        line = line.strip().split(sep=",")
        if line[0] == chromosome:
            if line[21] == band:
                pop = line[1]
                pop_ = "_" + pop
                start = int(line[10])
                end = int(line[11])
                xmin = min(xmin, start)
                xmax = max(xmax, end)
                clade_id = int(line[5])
                left_mut = int(line[14])
                right_mut = int(line[15])
                chunkindex = int(line[18])
                treeindex = int(line[19])
                node_id = int(line[20])
                print(node_id, pop)
                cladesize = int(line[8])
                p1 = float(line[6])
                p2 = float(line[7])

                ts = tszip.decompress(
                    trees_loc
                    + chromosome
                    + "/1000GP_Phase3_mask_prene_"
                    + chromosome
                    + "_chunk"
                    + str(chunkindex)
                    + ".trees.tsz"
                )

                node_found = False
                with open(time_file, "r") as tfile:
                    for tline in tfile:
                        tline = tline.strip().split(sep=",")
                        node_id_ = int(tline[1])
                        if node_id == node_id_:
                            node_found = True
                            mutations = [
                                m for m in tline[23].strip().split(sep=" ")
                            ]
                            if mutations[0] != "":
                                mutations = set(int(m) for m in mutations)
                                focal_mut = int(tline[12])
                                r2 = float(tline[5])
                                if r2 >= 0.95:
                                    chunkindex_ = int(tline[2])
                                    if chunkindex_ != chunkindex:
                                        chunkindex = chunkindex_
                                        ts = tszip.decompress(
                                            trees_loc
                                            + chromosome
                                            + "/1000GP_Phase3_mask_prene_"
                                            + chromosome
                                            + "_chunk"
                                            + str(chunkindex)
                                            + ".trees.tsz"
                                        )
                                    t = ts.at(focal_mut)
                                    if t.index == 1001:
                                        chunkindex += 1
                                        ts = tszip.decompress(
                                            trees_loc
                                            + chromosome
                                            + "/1000GP_Phase3_mask_prene_"
                                            + chromosome
                                            + "_chunk"
                                            + str(chunkindex)
                                            + ".trees.tsz"
                                        )
                                        t = ts.at(focal_mut)
                                    done_nodes = set()
                                    done_mutations = set()
                                    for mut in t.mutations():
                                        mp = int(ts.site(mut.site).position)
                                        if mp in mutations:
                                            done_mutations.add(mp)
                                            if mut.node not in done_nodes:
                                                counts = str(node_id) + "," + pop
                                                samples = [s for s in t.samples(mut.node)]
                                                pop_samples = [s for s in samples if md[s]["group"] == pop]
                                                for pp in populations:
                                                    freq_high = np.sum([1 for s in samples if md[s]["group"] == pp]) / num_samples[pp]
                                                    counts += "," + str(freq_high)
                                                samples_ = [s for s in t.samples(t.mrca(*pop_samples))]
                                                for pp in populations:
                                                    freq_low = np.sum([1 for s in samples_ if md[s]["group"] == pp]) / num_samples[pp]
                                                    counts += "," + str(freq_low)
                                                ffile.write(counts + "," + str(mp) + "\n")
                                                done_nodes.add(mut.node)
                                    if len(done_mutations) != len(mutations):
                                        print(t.index, t.interval)
                                        print(done_mutations)
                                        print(mutations)
                                        sys.exit("MUTATIONS NOT FOUND")
                if not node_found:
                    sys.exit("NODE NOT FOUND")

with open(info_file, "a") as ifile:
    ifile.write(band)
    print("Spans centromere:")
    if (xmin < centromeres[chromosome]) and (xmax > centromeres[chromosome]):
        print("yes")
        ifile.write(",yes")
    else:
        print("no")
        ifile.write(",no")
    xmin = xmin - 100000
    xmax = xmax + 100000

    dupstr_direct = ""
    dupstr_inverted = ""
    done_dups = []
    with open(superdups_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            chr1 = parts[1]
            chr2 = parts[7]
            if chr1 == chromosome or chr2 == chromosome:
                start1 = int(parts[2])
                end1 = int(parts[3])
                start2 = int(parts[8])
                end2 = int(parts[9])
                inv = parts[6]
                c1 = chromosome == chr1 and end1 >= xmin and start1 <= xmax
                c2 = chromosome == chr2 and end2 >= xmin and start2 <= xmax
                # print(start1, end1, start2, end2, chr1, chr2, chromosome, c1, c2)
                if c1 and (start1, end1) not in done_dups:
                    if c2 and (start2, end2) not in done_dups:
                        if inv == "-":
                            # inverted repeats
                            dupstr_inverted += str(start1) + "-" + str(end1) + "(" + str(end1-start1) + ")->"
                            dupstr_inverted += str(start2) + "-" + str(end2) + "(" + str(end2 - start2) + ");"
                        else:
                            dupstr_direct += str(start1) + "-" + str(end1) + "(" + str(end1-start1) + ")->"
                            dupstr_direct += str(start2) + "-" + str(end2) + "(" + str(end2 - start2) + ");"
                        done_dups.append((start2, end2))
                    done_dups.append((start1, end1))
    print("Seg dups:")
    ifile.write("," + dupstr_direct)
    ifile.write("," + dupstr_inverted)
    print(dupstr_direct)
    print(dupstr_inverted)

    cnv_start = math.inf
    cnv_end = -1
    with open(cnv_file, "r") as file:
        for line in file:
            line = line.split(sep=",")
            if line[0] == chromosome:
                p1 = int(line[3])
                p2 = int(line[4])
                if p1 < xmax and p2 > xmin:
                    cnv_start = min(cnv_start, p1)
                    cnv_end = max(cnv_end, p2)
    print("DoLoReS CNVs:")
    if cnv_end != -1:
        print(str(cnv_start) + "-" + str(cnv_end) + "(" + str(cnv_end - cnv_start) + ")")
        ifile.write("," + str(cnv_start) + "-" + str(cnv_end) + "(" + str(cnv_end - cnv_start) + ")")
    else:
        print("-")
        ifile.write(",")
    ifile.write("\n")
