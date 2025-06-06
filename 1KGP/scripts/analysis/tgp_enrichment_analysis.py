import sys
import os
import numpy as np
import random

from collections import defaultdict

sys.path.append("..")
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

# ch = int(sys.argv[1])

results_dir = "files/"
infile = results_dir + "SNP_samples_final.txt"

data = defaultdict(list)
observed = [0, 0, 0, 0, 0, 0, 0]
matched = {}
bands = set()

# Read in SNPs
with open(infile, "r") as file:
    next(file)
    for line in file:
        line = line.strip().split("\t")
        band = line[0]
        if band not in matched:
            matched[band] = [0, 0, 0, 0, 0, 0, 0, 0]
        bands.add(band)
        focal = int(line[1])
        pos = int(line[2])
        selection = int(int(line[11]) > 0)
        selection_50kb = int(int(line[12]) > 0)
        selection_100kb = int(int(line[13]) > 0)
        selection_200kb = int(int(line[14]) > 0)
        geneoverlaps = int(int(line[15]) > 0)
        meiotic = int(int(line[16]) > 0)
        meiotic_strong = int(int(line[17]) > 0)
        geneoverlaps_strong = int(int(line[18]) > 0)
        gwas = line[19]
        if gwas != "0":
            gwas = 1
        else:
            gwas = 0
        eqtls = line[20]
        if eqtls != "0":
            eqtls = eqtls.strip().split(");")
            eqtl_genes = []
            for eqtl in eqtls:
                eqtl = eqtl.split("(")
                gene = eqtl[0]
                eqtl_list = eqtl[1]
                eqtl_list = eqtl_list.split(";")
                for ee in eqtl_list:
                    if ee[0] == "e":
                        ee = ee.split(":")
                        coef = float(ee[1])
                        # print(band, gene, ee, coef)
                        if coef >= 0.5:
                            eqtl_genes.append((gene, coef))
            if len(eqtl_genes) > 0:
                eqtl = 1
            else:
                eqtl = 0
        else:
            eqtl = 0
        # print(band, focal, pos, selection, geneoverlaps, meiotic, meiotic_strong, geneoverlaps_strong, gwas, eqtl)
        data[band].append((selection, geneoverlaps, meiotic, meiotic_strong, geneoverlaps_strong, gwas, eqtl))
        if focal == 1:
            observed[0] += selection
            observed[1] += geneoverlaps
            observed[2] += meiotic
            observed[3] += meiotic_strong
            observed[4] += geneoverlaps_strong
            observed[5] += gwas
            observed[6] += eqtl
            if geneoverlaps_strong > 0:
                print(band)
        matched[band][0] += selection
        matched[band][1] += geneoverlaps
        matched[band][2] += meiotic
        matched[band][3] += meiotic_strong
        matched[band][4] += geneoverlaps_strong
        matched[band][5] += gwas
        matched[band][6] += eqtl
        matched[band][7] += 1

for band in matched:
    for i in range(7):
        matched[band][i] = matched[band][i]/matched[band][7]
print(matched)
print(observed)
print(len(matched))

p = [0, 0, 0, 0, 0, 0, 0, 0]
RR = 1000
for R in range(RR):
    results = [0, 0, 0, 0, 0, 0, 0]
    for band in bands:
        # print(data[band])
        selection, geneoverlaps, meiotic, meiotic_strong, geneoverlaps_strong, gwas, eqtl = random.sample(data[band], 1)[0]
        results[0] += selection
        results[1] += geneoverlaps
        results[2] += meiotic
        results[3] += meiotic_strong
        results[4] += geneoverlaps_strong
        results[5] += gwas
        results[6] += eqtl
    p[0] += int(observed[0] <= results[0])
    p[1] += int(observed[1] <= results[1])
    p[2] += int(observed[2] <= results[2])
    p[3] += int(observed[3] <= results[3])
    p[4] += int(observed[4] <= results[4])
    p[5] += int(observed[5] <= results[5])
    p[6] += int(observed[6] <= results[6])

labels = ["selection", "gene overlaps", "meiotic", "meiotic strong", "gene overlaps strong", "gwas", "eqtl"]
for i in range(len(labels)):
    print("- "*50)
    print(labels[i], p[i]/RR)

    P = [item[i] for item in matched.values()]
    print("Observed", observed[i])
    print("Expected", sum(P))
    pb = poibin.PoiBin(P)
    x = [j for j in range(len(matched) + 1)]
    pvals = pb.pval(x)
    cdfs = pb.cdf(x)
    for j in range(len(matched) + 1):
        if j == observed[i]:
            s = "--->"
        else:
            s = "    "
        print(s, x[j], pvals[j], cdfs[j])

