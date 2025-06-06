import sys
import os
import gzip
import glob
import numpy as np
import scipy as sp
from scipy.stats import pearsonr
import pickle
import tqdm
import tszip
import math

sys.path.append("..")
import dolores.simulations
import dolores.cladespans

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

chromosome = "chr8"
ch = 8
region = [6922488, 12573597]
focal_snps = [8633548]
populations = ["EUR", "AFR", "EAS", "SAS", "AMR"]

results_dir = "results_1000GP/results_varNe/tgp_" + chromosome + "/"
results_file = results_dir + "clades_pvalues_" + chromosome + "_EUR_HapMapII_GRCh37_events.csv"
trees_loc = "trees/relate_1000GP/trees/chunks/" + chromosome + "/"

# md, pops, groups = dolores.simulations.read_in_pop_metadata(
#     "trees/relate_1000GP/data/1000GP_Phase3.poplabels"
# )
# samp_to_ind_all = {}
# inds = {}
# i = 0
# for key, val in md.items():
#     if val["ID"] not in inds:
#         inds[val["ID"]] = i
#         i += 1
#     samp_to_ind_all[key] = inds[val["ID"]]


def get_bitset(n_set, samp_to_ind):
    e = [0] * int(len(samp_to_ind) / 2)
    for n in n_set:
        e[samp_to_ind[n]] += 1
    return e

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

for focal_snp in focal_snps:
    chunk_left = get_chunk(ch, region[0])
    chunk_right = get_chunk(ch, region[1])
    chunk = get_chunk(ch, focal_snp)
    print("Chunks:", chunk_left, chunk, chunk_right)

    for pop in populations:
        ts = tszip.decompress(trees_loc + "1000GP_Phase3_mask_prene_chr8_chunk" + str(chunk) + "_" + pop + ".trees.tsz")
        samp_to_ind = pickle.load(
            open(
                trees_loc
                + "1000GP_Phase3_mask_prene_"
                + chromosome
                + "_"
                + pop
                + "_samp_to_ind.pickle",
                "rb",
            )
        )

        for m in ts.mutations():
            p = ts.site(m.site).position
            if abs(p - focal_snp) < 2:
                print(pop, focal_snp, p, m)
                t = ts.at(p)
                inv_s = [s for s in t.samples(m.node)]
                non_s = [s for s in t.samples() if s not in inv_s]
                print(pop, "Frequency", len(inv_s)/num_samples[pop], len(non_s)/num_samples[pop])
                break

        bitset = get_bitset(set(inv_s), samp_to_ind)
        clades = dolores.cladedurations_.read_from_file("results_1000GP/results_varNe/tgp_" + chromosome + "/tgp_" + chromosome + "_" + pop + "_merged1.0_mutlimit0.0_HapMapII_GRCh37")

        chunkmin = chunk_left
        ts = tszip.decompress(trees_loc + "1000GP_Phase3_mask_prene_chr8_chunk" + str(chunkmin) + "_" + pop + ".trees.tsz")
        t = ts.at(region[0])
        all_nodes = []
        found_clades = []
        found_clades_r = []
        while region[0] <= t.interval[1] and t.interval[0] <= region[1]:
            if t.num_roots == 1:
                top_r = []
                top_nodes = []
                for n in t.nodes():
                    if n != t.root and not t.is_sample(n):
                        samples_ = {s for s in t.samples(n)}
                        b = get_bitset(samples_, samp_to_ind)
                        if len(b) != len(bitset):
                            print(len(b), np.sum(b))
                            sys.exit("len(b) != len(bitset)")
                        r = pearsonr(b, bitset)[0]
                        if r > 0.8:
                            top_r.append(r)
                            top_nodes.append(n)
                if len(top_r) > 0:
                    print(chunkmin, t.index, t.interval, max(top_r))
                all_nodes.extend(top_nodes)
                for i in range(len(top_nodes)):
                    n = top_nodes[i]
                    if n in clades.nodeid:
                        m = [i for i in range(clades.num) if clades.nodeid[i] == n and clades.chunkindex[i] == chunkmin and clades.treeindex[i] == t.index]
                        if len(m) > 1:
                            sys.exit("More than one clade found for node")
                        if len(m) == 1:
                            found_clades.append(m[0])
                            found_clades_r.append(top_r[i])
                if t.index == ts.num_trees - 2:
                    chunkmin += 1
                    ts =  tszip.decompress(trees_loc + "1000GP_Phase3_mask_prene_chr8_chunk" + str(chunkmin) + "_" + pop + ".trees.tsz")
                    t = ts.first()
            t.next()
        print(len(all_nodes))
        print(len(found_clades))

        with open("results_1000GP/" + chromosome + ".csv", "a") as file:
            for i in range(len(found_clades)):
                n = found_clades[i]
                r = found_clades_r[i]
                file.write(
                    pop
                    + ","
                    + str(n)
                    + ","
                    + str(clades.cladesize[n])
                    + ","
                    + str(r)
                    + ","
                    + str(clades.num_mutations[n])
                    + ","
                    + str(clades.start[n])
                    + ","
                    + str(clades.end[n])
                    + ","
                    + str(clades.duration[n])
                    + "\n"
                )








