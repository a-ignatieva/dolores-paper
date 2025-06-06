import sys
import os
import gzip
import glob
import stdpopsim
import numpy as np
import scipy as sp
from scipy.stats import pearsonr
import tqdm
import tszip
import math

sys.path.append("..")
import dolores.simulations

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

# ch = int(sys.argv[1])

results_dir = "results_1000GP/"
results_file = results_dir + "clades_pvalues_all_events_annotated_varNe.csv"
sel_file = "files/Selection_Summary_Statistics.tsv.gz"
trees_loc = "trees/relate_1000GP/trees/chunks/"
annot_file = "files/ncbiRefSeqCurated.txt"
gtex_file = "files/gtexGeneV8.txt.gz"
gtex_tissues_file = "files/gtexTissueV8.csv"
nodes_file = "files/bands_to_freq_nodes.csv"
meiosis_file = "files/meiotic_top_genes_by_gtex_exp.csv"
mut_summary_file = results_dir + "SNP_summary_"
outfile = results_dir + "SNP_samples_final.csv"

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
    "1q25.1-2",  # low freq in EUR
    "12q24.13",  # low freq in EUR
    "15q13.3",  # low freq in EUR
    "2p23.1-2",  # freq near 1 in EUR
    "2p14",  # freq near 1 in EUR
    "2q32.1",  # freq near 1 in EUR
    "7p22.2",  # freq near 1 in EUR
    "8p11.1",  # freq near 1 in EUR
    "14q23.3",  # freq near 1 in EUR
    "16p12.2",  # freq near 1 in EUR
    "6p11.1",  # pericentromeric
    "8p11.1",  # pericentromeric
    "10q11.1",  # pericentromeric
}

species = stdpopsim.get_species("HomSap")

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
populations = ["EUR", "AFR", "SAS", "EAS", "AMR"]
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

def get_last_chunk(ch):
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
                chunk = int(line[1])
    return chunk

# Nodes used for frequencies
nodes = {}
with open(nodes_file, "r") as file:
    for line in file:
        line = line.strip().split(",")
        nodes[line[0]] = int(line[1])
# print(nodes)

# Meiosis genes
meiosis = []
with open(meiosis_file, "r") as file:
    next(file)
    for line in file:
        line = line.strip().split(",")
        meiosis.append(line[0])
# print(meiosis)

# Region coords
ch_to_bands = {i:[] for i in range(1, 23)}
band_sizes = {}
with open(results_file, "r") as file:
    for l, line in enumerate(file):
        if l > 0:
            line = line.strip().split(",")
            chromosome = int(line[0][3:])
            band = line[21]
            if band not in ignore:
                if band not in ch_to_bands[chromosome]:
                    ch_to_bands[chromosome].append(band)
                left = int(line[10])
                right = int(line[11])
                if band in band_sizes:
                    band_sizes[band][0] = min(band_sizes[band][0], left)
                    band_sizes[band][1] = max(band_sizes[band][1], right)
                else:
                    band_sizes[band] = [left, right]

# Selection SNPs
selection_snps = {i: [] for i in range(23)}
with gzip.open(sel_file, "rt") as file:
    for l, line in enumerate(file):
        if line[0] != "#" and line[0] != "C":
            line = line.strip().split("\t")
            X = float(line[10])
            if abs(X) > 5.45:
                chromosome = int(line[0])
                pos = int(line[1])
                selection_snps[chromosome].append(pos)

# Gene coords
gene_regions = {band: 0 for band in band_sizes}
all_genes = {i: {} for i, l in chromosome_lengths.items()}
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
            if gene in all_genes[ch]:
                all_genes[ch][gene] = (min(start, all_genes[ch][gene][0]), max(end, all_genes[ch][gene][1]))
            else:
                all_genes[ch][gene] = (start, end)
        except ValueError:
            valueerrors += 1
            continue
# print("ValueError", valueerrors, lines, round(100*valueerrors/lines, 2))
count_genes = 0
for ch in all_genes:
    print(ch, len(all_genes[ch]))
    count_genes += len(all_genes[ch])
print(len(meiosis), count_genes)

#
# with open(outfile, "w") as ofile:
#     for band in nodes:
#         if band not in ignore:
#             print(band)
#             left = band_sizes[band][0]
#             right = band_sizes[band][1]
#             band_size = right - left + 1
#             if "p" in band:
#                 ch = int(band.strip().split("p")[0])
#             else:
#                 ch = int(band.strip().split("q")[0])
#             contig = species.get_contig(chromosome="chr" + str(ch), genetic_map="HapMapII_GRCh37")
#             rmap = contig.recombination_map
#             gensize1 = 100 * (rmap.get_cumulative_mass(right) - rmap.get_cumulative_mass(left))
#             snps_freq = []
#             snps = []
#             with open("results_1000GP/chr" + band + "_freq.csv", "r") as file:
#                 for line in file:
#                     line = line.strip().split(",")
#                     snps_freq.append([float(l) for l in line[2:7]])
#                     snps.append(int(line[-1]))
#             middle = (np.min(snps) + np.max(snps))/2  # middle of region
#             snps_dist = [abs(s - middle) for s in snps]
#             snp = [snps[i] for i in range(len(snps)) if snps_dist[i] == np.min(snps_dist)][0]
#             snp_freq = [snps_freq[i] for i in range(len(snps)) if snps_dist[i] == np.min(snps_dist)][0]
#             # print(snp, np.min(snps), np.max(snps), len(snps), snp_freq)
#
#             with open(mut_summary_file + "chr" + str(ch) + ".txt", "r") as file:
#                 for line in file:
#                     line = line.strip().split(",")
#                     pos = int(float(line[1]))
#                     if pos == snp:
#                         ld = [float(f) for f in line[7:12]]
#                         snp_ld = ld
#                         break
#
#             matched_snps = {snp: [snp_freq, snp_ld, gensize1, 0, [left, right]]}
#             with open(mut_summary_file + "chr" + str(ch) + ".txt", "r") as file:
#                 for line in file:
#                     line = line.strip().split(",")
#                     frequencies = [float(f) for f in line[2:7]]
#                     ld = [float(f) for f in line[7:12]]
#                     # diff_freq = [abs(f1 - f2) for f1, f2 in zip(snp_freq, frequencies)]  # won't use this - just want to match EUR freq
#                     diff_freq1 = abs(np.sum(frequencies) - np.sum(snp_freq))  # avg frequency
#                     diff_freq2 = abs(frequencies[0] - snp_freq[0])  # frequency in europeans
#                     diff_ld = abs(ld[-1] - snp_ld[-1])  # cM in region 10kb, 50kb, 100kb, 500kb, 1Mb nearby
#                     if diff_freq1 < 0.1 and diff_freq2 < 0.05 and diff_ld < 0.2:
#                         pos = int(float(line[1]))
#                         if np.min([abs(p - pos) for p in matched_snps]) > 500000:  # Need far enough away from focal SNP
#                             if pos == snp:
#                                 region = [left, right]
#                             else:
#                                 region = [int(pos - band_size / 2), int(pos + band_size / 2)]
#                             check_mask = False
#                             for m in mask[ch]:
#                                 if region[1] >= m[0] and region[0] <= m[1]:
#                                     check_mask = True
#                             if not check_mask:
#                                 gensize2 = 100 * (
#                                     rmap.get_cumulative_mass(region[1]) - rmap.get_cumulative_mass(region[0])
#                                 )
#                                 diff_gensize = abs(gensize1 - gensize2)
#                                 if diff_gensize < 0.2:
#                                     matched_snps[pos] = [frequencies, ld, gensize2, diff_gensize, region]
#             print(band, len(matched_snps))
#             diffs = np.array([vv[3] for vv in matched_snps.values()])
#             poss = [pos for pos in matched_snps.keys()]
#             if len(diffs) > 10:  # will choose 10 regions with closest genetic size to region
#                 idx = np.argpartition(diffs, 10)  # get indices of 10 smallest elements
#                 poss = [poss[i] for i in idx[0:10]]  # get corresponding positions
#             for m, (v1, v2, gensize2, _, region) in matched_snps.items():
#                 if m in poss:
#                     s = 0
#                     if m == snp:
#                         s = 1
#
#                     selection_snps_count = selection_snps_count_50kb = selection_snps_count_100kb = 0
#                     selection_snps_count_200kb = 0
#                     for pos in selection_snps[ch]:
#                         if region[0] <= pos <= region[1]:
#                             selection_snps_count += 1
#                         if region[0]-50000 <= pos <= region[1]+50000:
#                             selection_snps_count_50kb += 1
#                         if region[0]-100000 <= pos <= region[1]+100000:
#                             selection_snps_count_100kb += 1
#                         if region[0]-200000 <= pos <= region[1]+200000:
#                             selection_snps_count_200kb += 1
#
#                     gene_overlaps = []
#                     meiosis_count = 0
#                     meiosis_strong_count = 0
#                     for gene, (start, end) in all_genes[ch].items():
#                         if start < region[1] and end > region[0]:
#                             a1 = region[1] - max(region[0], start)
#                             a2 = region[1] - min(region[1], end)
#                             A = (a1 - a2) / (region[1] - region[0])
#                             B = (a1 - a2) / (end - start)
#                             gene_overlaps.append((gene, A, B))
#                             if gene in meiosis:
#                                 meiosis_count += 1
#                     gene_overlaps_strong = []
#                     for gene, A, B in gene_overlaps:
#                         if A > 0.8 and B > 0.8:
#                             gene_overlaps_strong.append((gene, A, B))
#                             if gene in meiosis:
#                                 meiosis_strong_count += 1
#                     ofile.write(
#                         band + ","
#                         + str(s) + ","
#                         + str(m) + ","
#                         + " ".join([str(round(vv, 2)) for vv in v1]) + ","
#                         + " ".join([str(round(vv, 2)) for vv in v2]) + ","
#                         + str(band_size) + ","
#                         + str(region[1] - region[0]) + ","
#                         + str(region[0]) + ","
#                         + str(region[1]) + ","
#                         + str(round(gensize1, 2)) + ","
#                         + str(round(gensize2, 2)) + ","
#                         + str(selection_snps_count) + ","
#                         + str(selection_snps_count_50kb) + ","
#                         + str(selection_snps_count_100kb) + ","
#                         + str(selection_snps_count_200kb) + ","
#                         + str(len(gene_overlaps)) + ","
#                         + str(meiosis_count) + ","
#                         + str(meiosis_strong_count) + ","
#                         + str(len(gene_overlaps_strong)) + "\n"
#                     )

# with open(mut_summary_file + "chr" + str(ch) + ".txt", "w") as file:
#     contig = species.get_contig(
#         chromosome="chr" + str(ch), genetic_map="HapMapII_GRCh37"
#     )
#     recm = contig.recombination_map
#     lastchunk = get_last_chunk(ch)
#     print(ch, lastchunk)
#     count = 0
#
#     with tqdm.tqdm(total = lastchunk + 1) as pbar:
#         for chunk in range(0, lastchunk + 1):
#             ts = tszip.decompress(
#                 trees_loc
#                 + "chr" + str(ch)
#                 + "/1000GP_Phase3_mask_prene_"
#                 + "chr" + str(ch)
#                 + "_chunk"
#                 + str(chunk)
#                 + ".trees.tsz"
#             )
#             for t in ts.trees():
#                 for m in t.mutations():
#                     samples = [s for s in t.samples(m.node)]
#                     if 10 <= len(samples) <= 5008-10:
#                         pos = ts.site(m.site).position
#                         if 1000000 < pos < recm.sequence_length - 1000000:
#                             summary = str(ch) + "," + str(pos)
#                             for pp in populations:
#                                 summary += "," + str(np.sum([1 for s in samples if md[s]["group"] == pp]) / num_samples[pp])
#                             for rsize in [10000, 50000, 100000, 500000, 1000000]:
#                                 gensize = recm.get_cumulative_mass(min(recm.sequence_length-1, pos + rsize/2)) - recm.get_cumulative_mass(max(0, pos - rsize/2))
#                                 summary += "," + str(gensize*100)
#                             summary += "\n"
#                             file.write(summary)
#                             count += 1
#             pbar.update(1)
#         print(count, "SNPs")

# band = "6p21.2"
# signif_snps = {}
# with open("results_1000GP/chr6p21.2_freq.csv", "r") as file:
#     for line in file:
#         line = line.strip().split(",")
#         signif_snps[int(line[-1])] = int(line[0])
# ldsnps = []

# gtex_tissues = []
# with open(gtex_tissues_file, "r") as file:
#     next(file)
#     for line in file:
#         line = line.strip().split(",")
#         gtex_tissues.append(line[1][1:-1])
# print(gtex_tissues)
# print(len(gtex_tissues))
#
# gtex = {}
# with gzip.open(gtex_file, "rt") as file:
#     for line in file:
#         line = line.strip().split()
#         scores = line[9].strip()[:-1].split(",")
#         scores = np.array([float(g) for g in scores])
#         top5 = np.argpartition(scores, -5)[-5:]
#         gtex[line[3]] = (line[7], [gtex_tissues[i] for i in top5], scores[top5])
# print(gtex["COMMD10"])
# print(len(gtex["COMMD10"][1]))

# chromosomes = {i:[] for i in range(1, 23)}
# bands = {}
# with open(results_file, "r") as file:
#     for l, line in enumerate(file):
#         if l > 0:
#             line = line.strip().split(",")
#             chromosome = int(line[0][3:])
#             band = line[21]
#             if band not in ignore:
#                 if band not in chromosomes[chromosome]:
#                     chromosomes[chromosome].append(band)
#                 left = int(line[10])
#                 right = int(line[11])
#                 if band in bands:
#                     bands[band][0] = min(bands[band][0], left)
#                     bands[band][1] = max(bands[band][1], right)
#                 else:
#                     bands[band] = [left, right]
# lengths = {band: bands[band][1] - bands[band][0] for band in bands}
#
# print("="*100)
# print("Regions:", len(bands))
# print("- "*50)


#
# for chromosome in range(1, 23):
#     print("="*100)
#     print("Chromosome", chromosome)
#     if len(chromosomes[chromosome]) != 0:
#         print("Selection SNPs:")
#         print(len(selection_snps[chromosome]))
#
#         print("Proportion of chromosome that is genes")
#         print(chromosome, np.sum([g[1] - g[0] for g in all_genes[chromosome].values()])/chromosome_lengths[chromosome])
#
#         contig = species.get_contig(chromosome="chr" + str(chromosome), genetic_map="HapMapII_GRCh37")
#         recombination_map = contig.recombination_map
#
#         for band in chromosomes[chromosome]:
#             print("- "*50)
#             print("Region:", band)
#             print("Genomic coords:", bands[band])
#             left = bands[band][0]
#             right = bands[band][1]
#             size = right - left
#             print("region length (bp)", right - left + 1)
#             gensize = recombination_map.get_cumulative_mass(right) - recombination_map.get_cumulative_mass(left)
#             print("region length (cM)", gensize * 100)
#
#             selection_snps_count = 0
#             for pos in selection_snps[chromosome]:
#                 if left <= pos <= right:
#                     selection_snps_count += 1
#             print("selection SNPs:", selection_snps_count)
#
#             print("computing gene overlaps...")
#             gene_overlaps = []
#             for gene, (start, end) in all_genes[chromosome].items():
#                 if start < right and end > left:
#                     a1 = right - max(left, start)
#                     a2 = right - min(right, end)
#                     A = (a1 - a2) / (right - left)
#                     B = (a1 - a2) / (end - start)
#                     gene_overlaps.append((gene, A, B))
#             for gene, A, B in gene_overlaps:
#                 s = "     "
#                 if A > 0.8 and B > 0.8:
#                     s = "---> "
#                 print(s, gene, round(A, 3), round(B, 3))


