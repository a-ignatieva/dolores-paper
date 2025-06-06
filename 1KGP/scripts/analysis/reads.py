#!/usr/bin/env python
# coding: utf-8

import os
import sys

import pysam
import numpy as np
import math
import random
import pickle
import tszip

def read_in_pop_metadata(ff):
    """
    Function to read in file with 4 columns (ID, POP, GROUP, SEX) and return a dictionary with the data,
    and group and population labels
    :param file:
    :return:
    """
    metadata = {}
    grps = set()
    pops = set()
    with open(ff, "r") as f:
        next(f)
        i = 0
        for line in f:
            line = line.strip().split()
            metadata[i] = {
                "ID": line[0],
                "population": line[1],
                "group": line[2],
                "sex": line[3],
            }
            metadata[i + 1] = {
                "ID": line[0],
                "population": line[1],
                "group": line[2],
                "sex": line[3],
            }
            if line[2] not in grps:
                grps.add(line[2])
            if line[1] not in pops:
                pops.add(line[1])
            i += 2
    return metadata, pops, grps

def run_reads(todo, done, indexer, ch, start, end):
    count_ = 0
    with (open(band + "_insert-sizes.txt", "a") as file_insert,
          open(band + "_coverage.txt", "a") as file_cov,
          open(band + "_done.txt", "a") as file_done):
        for p, (u, pop, pop_) in indexer.items():
            if p in todo and p not in done:
                print(count_, count, p)

                file_insert.write(
                    p
                    + ","
                    + pop
                    + ","
                    + pop_
                    + ","
                    + "NA"
                    + ","
                    + "NA"
                    + ","
                    + "NA"
                    + "\n"
                )

                os.system(
                    "samtools view -h -b "
                    + u
                    + " "
                    + str(ch)
                    + ":"
                    + str(start)
                    + "-"
                    + str(end)
                    + " > " + band + "_test.bam"
                )
                os.system("samtools index " + band + "_test.bam")

                skip = False
                try:
                    samfile = pysam.Samfile(band + "_test.bam", 'rb')
                except:
                    skip = True
                if not skip:
                    for read in samfile.fetch():
                        if read.is_mapped and read.mate_is_mapped:
                            if read.is_paired:
                                if read.reference_id + 1 == ch and read.next_reference_id + 1 == ch:
                                    discordant = 0
                                    if (read.is_reverse and read.mate_is_reverse) or (
                                            not read.is_reverse and not read.mate_is_reverse):
                                        discordant = 1
                                    if discordant == 1 or (1000 <= abs(read.next_reference_start - read.reference_start) <= (end - start) + 1000000):
                                        file_insert.write(
                                            p
                                            + ","
                                            + pop
                                            + ","
                                            + pop_
                                            + ","
                                            + str(read.reference_start)
                                            + ","
                                            + str(read.next_reference_start)
                                            + ","
                                            + str(abs(read.next_reference_start - read.reference_start))
                                            + ","
                                            + str(discordant)
                                            + "\n"
                                        )

                    file_cov.write(
                        p
                        + ","
                        + pop
                        + ","
                        + pop_
                    )
                    counts = np.zeros(len(positions))
                    for i, pileupcolumn in enumerate(samfile.pileup(str(ch), start, end)):
                        if start <= pileupcolumn.pos < end:
                            counts[positions[pileupcolumn.pos]] += pileupcolumn.n
                    for i in range(int(len(counts) / 1000)):
                        ind1 = 1000 * i
                        ind2 = min(len(counts) - 1, 1000 * (i + 1))
                        ct = np.sum(counts[ind1:ind2]) / 1000
                        file_cov.write("," + str(ct))
                    file_cov.write("\n")
                    file_done.write(p + "\n")

                    samfile.close()
                os.system("rm " + p + "*.bai")

            count_ += 2

populations = {
    "CEU": "EUR", "TSI": "EUR", "GBR": "EUR", "FIN": "EUR", "IBS": "EUR",
    'MSL': "AFR", 'ESN': "AFR", 'ACB': "AFR", 'YRI': "AFR", 'LWK': "AFR", 'ASW': "AFR", 'GWD': "AFR",
    'KHV': "EAS", 'CDX': "EAS", 'JPT': "EAS", 'CHS': "EAS", 'CHB': "EAS",
    'ITU': "SAS", 'BEB': "SAS", 'GIH': "SAS", 'STU': "SAS", 'PJL': "SAS",
    'PEL': "AMR", 'CLM': "AMR", 'MXL': "AMR", 'PUR': "AMR",
}
results_file = "files/clades_pvalues_all_events_annotated_varNe.csv"
trees_loc = "trees/relate_1000GP/trees/chunks/"
reads_loc = "reads/"
# nsamples = 50

band = sys.argv[1]
if "p" in band:
    ch = int(band.split("p")[0])
else:
    ch = int(band.split("q")[0])
chromosome = "chr" + str(ch)
print(ch, band)

md, pops, groups = read_in_pop_metadata(
    "files/1000GP_Phase3.poplabels"
)
print(pops, groups)

start = math.inf
end = -1

with open(results_file, "r") as file:
    for line in file:
        line = line.strip().split(",")
        if line[21] == band:
            start = min(start, int(line[10]))
            end = max(end, int(line[11]))
start = int(start - 100000)
end = int(end + 100000)
positions = {s: i for i, s in enumerate(range(start, end))}
print(start, end)
with open(band + "_info.txt", "w") as file:
    file.write(str(start) + "," + str(end) + "\n")

# Collect the URLs
print("Collecting URLs...")
populations_counts = {"EUR": 0, "AFR": 0, "SAS": 0, "EAS": 0, "AMR": 0}
indexer = {}
count = 0
with open(reads_loc + "igsr_Phase 3.tsv", "r") as index_file:
    for l, line in enumerate(index_file):
        if l > 0:
            line = line.strip().split(sep="\t")
            u = line[0]
            check1 = u.split(sep="/")
            if check1[8] == "alignment":
                check2 = check1[9].split(sep=".")
                if check2[1] == "mapped" and check2[-1] != "bai":
                    # print(u, check2[0], check2[4])
                    pop = check2[4]
                    pop_ = populations[pop]
                    indexer[check2[0]] = (u, pop, pop_)
                    populations_counts[pop_] += 2
                    count += 2
print(count, populations_counts)

done = set()
if os.path.exists(band + "_done.txt"):
    with open(band + "_done.txt", "r") as file:
        for line in file:
            line = line.strip()
            done.add(line)

todo = [v['ID'] for v in md.values()]
run_reads(todo, done, indexer, ch, start, end)

# with open(results_file, "r") as file:
#     for line in file:
#         line = line.strip().split(",")
#         if line[21] == band:
#             pop = line[1]
#             chunkindex = int(line[18])
#             treeindex = int(line[19])
#             node_id = int(line[20])
#             print(node_id)
#             samp_to_ind = pickle.load(
#                 open(
#                     trees_loc
#                     + chromosome
#                     + "/1000GP_Phase3_mask_prene_"
#                     + chromosome
#                     + "_"
#                     + pop
#                     + "_samp_to_ind.pickle",
#                     "rb",
#                 )
#             )
#             ts = tszip.decompress(
#                 trees_loc
#                 + chromosome
#                 + "/1000GP_Phase3_mask_prene_"
#                 + chromosome
#                 + "_chunk"
#                 + str(chunkindex)
#                 + "_"
#                 + pop
#                 + ".trees.tsz"
#             )
#             mapp = pickle.load(
#                 open(
#                     trees_loc
#                     + chromosome
#                     + "/1000GP_Phase3_mask_prene_"
#                     + chromosome
#                     + "_"
#                     + pop
#                     + "_map.pickle",
#                     "rb",
#                 )
#             )
#
#             t = ts.at_index(treeindex)
#             samples = {s for s in t.samples(node_id)}
#             print("Number of samples in clade:", len(samples))
#             ind_index = {}
#             for v in md.values():
#                 ind_index[v['ID']] = 0
#             for s in samples:
#                 ss = np.where(mapp == s)[0][0]
#                 ind = md[ss]['ID']
#                 ind_index[ind] += 1
#
#             ind_list = [i for i in ind_list if ind_index[i] == 0]
#             hom = [i for i in ind_index if ind_index[i] == 2]
#             print("Number of homozygotes:", len(hom))
#             todo = random.sample(hom, min(len(hom), 50))
#             print("Sample:", todo)
#
#             run_reads(todo, indexer, node_id, ch, start, end)

# todo = random.sample(ind_list, min(len(ind_list), 50))
# run_reads(todo, indexer, -1, ch, start, end)

os.system("rm " + band + "_test.*")
