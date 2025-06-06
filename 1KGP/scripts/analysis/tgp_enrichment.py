import sys
import os
import gzip
import glob
import numpy as np
import scipy as sp
from scipy.stats import pearsonr
import tqdm
import tszip
import math

sys.path.append("..")
import dolores.simulations
import poibin

chromosome_lengths = {
    1:249250621,
    2:243199373,
    3:198022430,
    4:191154276,
    5:180915260,
    6:171115067,
    7:159138663,
    8:146364022,
    9:141213431,
    10:135534747,
    11:135006516,
    12:133851895,
    13:115169878,
    14:107349540,
    15:102531392,
    16:90354753,
    17:81195210,
    18:78077248,
    19:59128983,
    20:63025520,
    21:48129895,
    22:51304566,
}
mask = {
    1:[(121484388, 142540862)],
    2:[(89622811, 89830437), (90536653, 91624982), (92324287, 95326256)],
    3:[(90504828, 93504880)],
    4:[(49660116, 52660482)],
    5:[(46405622, 49405642)],
    6:[(58780017, 61880193)],
    7:[(58054317, 61054338)],
    8:[(43838850, 46846021)],
    9:[(42613299, 42827867), (47266441, 65471163), (70217995, 70478235)],
    10:[(39154935, 42354939)],
    11:[(50783762, 51092947), (51593465, 54695442)],
    12:[(34856677, 37856749)],
    13:[(0, 19020067)],
    14:[(0, 19000019)],
    15:[(0, 20000027), (21398720, 21885385)],
    16:[(35285789, 46410209)],
    17:[(22262824, 25263011)],
    18:[(15410715, 18510899)],
    19:[(24629288, 27731820)],
    20:[(26319554, 29419679)],
    21:[(0, 9411210), (11187995, 14338149)],
    22:[(0, 16050343)],
}

results_dir = "files/"
results_file = results_dir + "clades_pvalues_all_events_annotated_varNe.csv"
sel_file = "files/Selection_Summary_Statistics.tsv.gz"
trees_loc = "trees/relate_1000GP/trees/chunks/"
annot_file = "files/ncbiRefSeqCurated.txt"

overlap1 = overlap2 = 0.8

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

ignore = {
    # "1q25.1-2",  # low freq in EUR
    # "12q24.13",  # low freq in EUR
    # "15q13.3",  # low freq in EUR
    "6p11.1",  # pericentromeric
    "8p11.1",  # pericentromeric
    "10q11.1",  # pericentromeric
}

def get_chunk(ch, pos):
    chunk = None
    with open(
        "trees/relate_1000GP/trees/chunks/treeinfo_chr"
        + str(ch)
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


chromosomes = {i:[] for i in range(1, 23)}
bands = {}
with open(results_file, "r") as file:
    for l, line in enumerate(file):
        if l > 0:
            line = line.strip().split(",")
            chromosome = int(line[0][3:])
            band = line[21]
            if band not in ignore:
                if band not in chromosomes[chromosome]:
                    chromosomes[chromosome].append(band)
                left = int(line[10])
                right = int(line[11])
                if band in bands:
                    bands[band][0] = min(bands[band][0], left)
                    bands[band][1] = max(bands[band][1], right)
                else:
                    bands[band] = [left, right]
lengths = {band: bands[band][1] - bands[band][0] for band in bands}
print("="*100)
print(len(bands))
print(chromosomes)
print("Region sizes")
print(lengths)
print("="*100)

#======================================================================================================
# SELECTION
#======================================================================================================

# for padding in [0, 200000, 500000, 1000000]:
#     print("- "*50)
#     print("Padding", padding)
#     interesting_bands = set()
#     sel_regions = {band: 0 for band in bands}
#     prevch = -1
#     with gzip.open(sel_file, "rt") as file:
#         for l, line in enumerate(file):
#             if line[0] != "#" and line[0] != "C":
#                 line = line.strip().split("\t")
#                 X = float(line[10])
#                 if abs(X) > 5.45:
#                     chromosome = int(line[0])
#                     pos = int(line[1])
#                     if prevch != chromosome:
#                         prevpos = -1
#                         prevch = chromosome
#                     for band in chromosomes[chromosome]:
#                         if bands[band][0] - padding <= pos <= bands[band][1] + padding:
#                             interesting_bands.add(band)
#                         L = lengths[band] + 2*padding
#                         if prevpos <= pos - L:
#                             sel_regions[band] += L
#                         else:
#                             sel_regions[band] += pos - prevpos
#                     prevpos = pos
#
#     total_regions = {band: 0 for band in bands}
#     for chromosome, length in chromosome_lengths.items():
#         for band in bands:
#             prevpos = 0
#             L = lengths[band]
#             for m in mask[chromosome]:
#                 total_regions[band] += m[0] - prevpos - L + 1
#                 prevpos = m[1]
#             total_regions[band] += length - prevpos - L + 1
#
#     P = []
#     for b in sel_regions.keys():
#         p = sel_regions[b]/total_regions[b]
#         # print(b, sel_regions[b], total_regions[b], round(p, 4))
#         P.append(p)
#
#     print(interesting_bands)
#     print("Observed", len(interesting_bands))
#     print("Expectation", sum(P))
#     pb = poibin.PoiBin(P)  # Poisson Binomial
#     x = [i for i in range(len(bands) + 1)]
#     pvals = pb.pval(x)
#     cdfs = pb.cdf(x)
#     for i in range(15):
#         if i == len(interesting_bands):
#             s = "--->"
#         else:
#             s = "    "
#         print(s, x[i], pvals[i], cdfs[i])

#======================================================================================================
# GENES
#======================================================================================================

for padding in [0]:
    print("- "*50)
    print("Padding", padding)
    interesting_bands = set()
    gene_regions = {band: 0 for band in bands}
    all_genes = {i: np.zeros(l) for i, l in chromosome_lengths.items()}
    valueerrors = 0
    lines = 0

    with open(annot_file, "r") as file:
        for line in file:
            lines += 1
            line = line.split()
            try:
                ch = int(line[2][3:])
                gene = line[12]
                start = int(line[4])
                end = int(line[5])
                all_genes[ch][start:(end+1)] = 1
            except ValueError:
                valueerrors += 1
                # print(line)
                # print(line[2])
                continue
    for ch, masks in mask.items():
        for m in masks:
            all_genes[ch][m[0]:m[1]] = 0
    print("ValueError", valueerrors, lines, round(100*valueerrors/lines, 2))

    for i, l in chromosome_lengths.items():
        print("Proportion of chromosome that is genes")
        print(i, l, np.sum(all_genes[i])/l)
        all_genes_summary = []

        for pos in range(1, l):
            if all_genes[i][pos] == 1:
                if all_genes[i][pos-1] == 0:
                    all_genes_summary.append([pos, pos])
                else:
                    all_genes_summary[-1][1] = pos
        minlength = math.inf
        maxlength = avglength = count = 0
        for start, end in all_genes_summary:
            minlength = min(minlength, end - start + 1)
            maxlength = max(maxlength, end - start + 1)
            avglength += end - start + 1
            count += 1
        print("min length of gene", minlength, "max length of gene", maxlength, "mean length of gene", avglength/count)

        for band in chromosomes[i]:
            print("")
            print(band)
            # print("Calculating overlaps...")
            # if np.sum(all_genes[i][(bands[band][0] - padding):(bands[band][1] + padding)]) > 0:
            #     interesting_bands.add(band)

            xmin = bands[band][0]
            xmax = bands[band][1]
            L = lengths[band] + 2*padding
            print("Calculating gene_regions...")
            prevend = -1
            prevsize = -1
            for start, end in all_genes_summary:
                K = end - start
                A = (min(xmax, end) - max(xmin, start)) / (end - start)
                B = (min(xmax, end) - max(xmin, start)) / (xmax - xmin)
                if A > overlap1 and B > overlap2:
                    interesting_bands.add(band)
                    print(band, A, B, start, end, xmin, xmax)
                if L >= (1 + end - start)*overlap1:
                    gene_regions[band] += max(0, L + 1 + K - 2 * max(K*overlap1, L * overlap2))
                # gene_regions[band] += L + 1 + end - start
                # if prevend != -1:
                #     gene_regions[band] -= max(0, prevend + L - start + 1 - overlap1*(prevsize + K))
                prevend = end
                prevsize = K

    print("Calculating total regions...")
    total_regions = {band: 0 for band in bands}
    for chromosome, length in chromosome_lengths.items():
        for band in chromosomes[chromosome]:
            prevpos = 0
            L = lengths[band] + 2*padding
            for m in mask[chromosome]:
                total_regions[band] += m[0] - prevpos - L + 1
                prevpos = m[1]
            total_regions[band] += length - prevpos - L + 1

    print("Calculating probabilities...")
    P = []
    for b in gene_regions.keys():
        p = min(1 - 1e-10, gene_regions[b]/total_regions[b])
        print(b, gene_regions[b], total_regions[b], gene_regions[b] > total_regions[b], round(p, 4))
        P.append(p)

    print("Interesting bands:")
    print(interesting_bands)
    print("Non-interesting bands:")
    print([band for band in bands if band not in interesting_bands])
    print("Observed", len(interesting_bands))
    print("Expectation", sum(P))
    pb = poibin.PoiBin(P)  # Poisson Binomial
    x = [i for i in range(len(bands) + 1)]
    pvals = pb.pval(x)
    cdfs = pb.cdf(x)
    for i in range(len(bands)):
        if i == len(interesting_bands):
            s = "--->"
        else:
            s = "    "
        print(s, x[i], pvals[i], cdfs[i])


