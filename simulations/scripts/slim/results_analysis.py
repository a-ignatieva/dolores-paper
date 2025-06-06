import os
import sys
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from sklearn.metrics import adjusted_rand_score

loc = "/Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/"

col_green = "#228833"
col_red = "#EE6677"
col_purp = "#AA3377"
col_blue = "#66CCEE"
col_yellow = "#CCBB44"
col_indigo = "#4477AA"
col_grey = "#BBBBBB"

dolores_c2_all = []
dolores_rand_all = []
dolores_overlap_all = []
invclust_c2_all_true = []
invclust_c2_all_large = []
asaph_rand_all = []
colors = [col_red, col_blue, col_green]

for il, invlength in enumerate([50000, 100000, 200000]):
    dolores_c2_results = []
    dolores_rand_results = []
    dolores_overlap_results = []
    invclust_c2_results_true = []
    invclust_c2_results_large = []
    asaph_rand_results = []
    dolores_fail = invclust_fail = 0
    breakpoints = [2500000 - invlength/2, 2500000 + invlength/2]
    ncount = 0
    for i in range(100):
        loc_ = loc + "slim_" + str(invlength) + "_50_" + str(i) + "/"
        go1 = go2 = go3 = True
        if ncount == 100:
            break
        if not os.path.exists(loc_ + "dolores_breakpoints.txt") or not os.path.exists(loc_ + "dolores_prediction.txt"):
            go1 = False
            dolores_fail += 1
            # print("dolores", invlength, i)
        else:
            with open(loc_ + "dolores_breakpoints.txt", "r") as file:
                for line in file:
                    dolores_breakpoints = [float(l) for l in line.strip().split(",")]
            dolores_carriers = np.zeros(50)
            with open(loc_ + "dolores_prediction.txt", "r") as file:
                for k, line in enumerate(file):
                    line = line.strip().split(",")
                    dolores_carriers[k] += float(line[1])
            dolores_clusters = np.zeros(50)
            with open(loc_ + "dolores_prediction.txt", "r") as file:
                for k, line in enumerate(file):
                    line = line.strip().split(",")
                    dolores_clusters[k] = int(float(line[1]) * 2)
            carriers = np.zeros(50)
            with open(loc_ + "carriers.txt", "r") as file:
                for k, line in enumerate(file):
                    line = line.strip().split(",")
                    carriers[k] += float(line[1])
            carrier_clusters = np.zeros(50)
            with open(loc_ + "carriers.txt", "r") as file:
                for k, line in enumerate(file):
                    line = line.strip().split(",")
                    carrier_clusters[k] = int(float(line[1])*2)
            # print(dolores_breakpoints)
            # print(dolores_carriers)
            # print(carriers)
            dolores_c2 = (stats.pearsonr(dolores_carriers, carriers)[0])**2
            dolores_rand = adjusted_rand_score(carrier_clusters, dolores_clusters)
            if dolores_c2 < 0.95 and dolores_rand > 0.95:
                print(invlength, i)
            # print(dolores_c2)
            dolores_overlap = max(0, min(breakpoints[1], dolores_breakpoints[2]) - max(breakpoints[0], dolores_breakpoints[1]))/invlength
            # print(dolores_overlap)
        if not os.path.exists(loc_ + "invclust.txt"):
            go2 = False
            invclust_fail += 1
            # print("invclust", invlength, i)
        else:
            invclust_c2 = np.zeros(4)
            with open(loc_ + "invclust.txt", "r") as file:
                next(file)
                for v, line in enumerate(file):
                    line = line.strip().split()
                    invclust_c2[v] = float(line[1])
            # print(invclust_c2)
        if os.path.exists(loc_ + "asaph_score.txt"):
            with open(loc_ + "asaph_score.txt", "r") as file:
                for line in file:
                    line = line.strip().split()
                    asaph_rand = float(line[0])
                    break
        else:
            go3 = False
        if go1:
            if go2 and go3:
                dolores_c2_results.append(dolores_c2)
                dolores_rand_results.append(dolores_rand)
                dolores_overlap_results.append(dolores_overlap)
                invclust_c2_results_true.append(invclust_c2[0])
                invclust_c2_results_large.append(invclust_c2[1])
                asaph_rand_results.append(asaph_rand)
                ncount += 1
            else:
                dolores_overlap_results.append(dolores_overlap)
    dolores_c2_all.append(dolores_c2_results)
    dolores_rand_all.append(dolores_rand_results)
    dolores_overlap_all.append(dolores_overlap_results)
    invclust_c2_all_true.append(invclust_c2_results_true)
    invclust_c2_all_large.append(invclust_c2_results_large)
    asaph_rand_all.append(asaph_rand_results)

    print(invlength, "dolores c2", np.mean(dolores_c2_results))
    print(invlength, "dolores rand", np.mean(dolores_rand_results))
    print(invlength, "dolores overlap", np.mean(dolores_overlap_results))
    print(invlength, "invclust c2", np.mean(invclust_c2_results_true))
    print(invlength, "invclust c2 large", np.mean(invclust_c2_results_large))
    print(invlength, "asaph rand", np.mean(asaph_rand_results))
    print(invlength, (len([w for w in range(len(dolores_c2_results)) if dolores_c2_results[w] >= invclust_c2_results_true[w]]))/len(dolores_c2_results))
    print(len(dolores_c2_results)/100)

plt.figure(figsize=(3.5, 3))
counts1 = []
counts2 = []
counts3 = []
counts4 = []
for i in range(3):
    plt.scatter(invclust_c2_all_true[i], dolores_c2_all[i], alpha = 0.5, color=colors[i])
    L = len(dolores_c2_all[i])
    count1 = set([j for j in range(L) if dolores_c2_all[i][j] >= 0.8])
    count2 = set([j for j in range(L) if invclust_c2_all_true[i][j] >= 0.8])
    counts1.append(len(count1 & count2) / L)
    counts2.append(len(count1 - count2) / L)
    counts3.append(len(count2 - count1) / L)
plt.text(x=0.83, y=0.85, s=round(np.mean(counts1), 2))
plt.text(x=0.65, y=0.85, s=round(np.mean(counts2), 2))
plt.text(x=0.83, y=0.7, s=round(np.mean(counts3), 2))
plt.text(x=0.65, y=0.7, s=round(1 - np.mean(counts1) - np.mean(counts2) - np.mean(counts3), 2))
plt.xlabel("invClust R^2 with ground truth")
plt.ylabel("DoLoReS R^2 with ground truth")
plt.hlines(xmin=-1, xmax=2, y=0.8, linewidth=0.5, color="black")
plt.vlines(ymin=-1, ymax=2, x=0.8, linewidth=0.5, color="black")
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
plt.savefig("results_genotypes_true.pdf", dpi=300, bbox_inches="tight")

plt.figure(figsize=(3.5, 3))
counts1 = []
counts2 = []
counts3 = []
counts4 = []
for i in range(3):
    plt.scatter(invclust_c2_all_large[i], dolores_c2_all[i], alpha = 0.5, color=colors[i])
    L = len(dolores_c2_all[i])
    count1 = set([j for j in range(L) if dolores_c2_all[i][j] >= 0.8])
    count2 = set([j for j in range(L) if invclust_c2_all_large[i][j] >= 0.8])
    counts1.append(len(count1 & count2) / L)
    counts2.append(len(count1 - count2) / L)
    counts3.append(len(count2 - count1) / L)
plt.text(x=0.83, y=0.85, s=round(np.mean(counts1), 2))
plt.text(x=0.65, y=0.85, s=round(np.mean(counts2), 2))
plt.text(x=0.83, y=0.7, s=round(np.mean(counts3), 2))
plt.text(x=0.65, y=0.7, s=round(1 - np.mean(counts1) - np.mean(counts2) - np.mean(counts3), 2))
plt.xlabel("invClust R^2 with ground truth")
plt.ylabel("DoLoReS R^2 with ground truth")
plt.hlines(xmin=-1, xmax=2, y=0.8, linewidth=0.5, color="black")
plt.vlines(ymin=-1, ymax=2, x=0.8, linewidth=0.5, color="black")
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
plt.savefig("results_genotypes_large.pdf", dpi=300, bbox_inches="tight")

plt.figure(figsize=(3.5, 3))
counts1 = []
counts2 = []
counts3 = []
counts4 = []
for i in range(3):
    plt.scatter(asaph_rand_all[i], dolores_rand_all[i], alpha = 0.5, color=colors[i])
    L = len(dolores_rand_all[i])
    count1 = set([j for j in range(L) if dolores_rand_all[i][j] >= 0.8])
    count2 = set([j for j in range(L) if asaph_rand_all[i][j] >= 0.8])
    counts1.append(len(count1 & count2) / L)
    counts2.append(len(count1 - count2) / L)
    counts3.append(len(count2 - count1) / L)
plt.text(x=0.83, y=0.85, s=round(np.mean(counts1), 2))
plt.text(x=0.65, y=0.85, s=round(np.mean(counts2), 2))
plt.text(x=0.83, y=0.7, s=round(np.mean(counts3), 2))
plt.text(x=0.65, y=0.7, s=round(1 - np.mean(counts1) - np.mean(counts2) - np.mean(counts3), 2))
plt.xlabel("asaph adj Rand with ground truth")
plt.ylabel("DoLoReS adj Rand with ground truth")
plt.hlines(xmin=-1, xmax=2, y=0.8, linewidth=0.5, color="black")
plt.vlines(ymin=-1, ymax=2, x=0.8, linewidth=0.5, color="black")
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
plt.savefig("results_genotypes_asaph.pdf", dpi=300, bbox_inches="tight")

plt.figure(figsize=(3.5, 3))
counts = []
for i in range(3):
    plt.hist(dolores_overlap_all[i], alpha = 0.5, color=colors[i], density=True, bins=10)
    counts.append(len([p for p in dolores_overlap_all[i] if p > 0.5])/len(dolores_overlap_all[i]))
plt.vlines(ymin=0, ymax=6, x=0.5, linewidth=0.5, color="black")
plt.text(x=0.51, y=5.5, s=round(np.mean(counts), 2))
plt.text(x=0.35, y=5.5, s=round(1 - np.mean(counts), 2))
plt.ylim(0, 6)
plt.xlabel("DoLoReS predicted region overlap")
plt.ylabel("Density")
plt.savefig("results_overlap.pdf", dpi=300, bbox_inches="tight")














