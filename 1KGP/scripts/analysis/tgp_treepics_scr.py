#!/usr/bin/env python
# coding: utf-8

import sys
import pickle

import tszip

import numpy as np
from cairosvg import svg2pdf

sys.path.append("..")
import argbd.simulations
import argbd.viz

col_green = "#228833"
col_red = "#EE6677"
col_purp = "#AA3377"
col_blue = "#66CCEE"
col_yellow = "#CCBB44"
col_indigo = "#4477AA"
col_grey = "#BBBBBB"

colours = {
    "AMR": col_blue,
    "AFR": "orange",
    "EAS": col_green,
    "EUR": col_indigo,
    "SAS": col_purp,
}

chromosome = sys.argv[1]
band = sys.argv[2]
positions = sys.argv[3]
node_id_list = sys.argv[4]
node_to_calculate = sys.argv[5]
if node_to_calculate != "None":
    node_to_calculate = int(node_to_calculate)
ind_to_highlight = sys.argv[6]

name = "results_1000GP/" + chromosome + band + "_tree_check_"
positions = [int(i) for i in positions.split(sep=",")]
node_id_list = [int(i) for i in node_id_list.split(sep=",")]

results_file = (
    "files/clades_pvalues_all_events_annotated_lowth_varNe.csv"
)
tree_index_list = {}
chunk_index_list = {}
populations_list = {}
with open(results_file, "r") as file:
    for line in file:
        line = line.strip().split(sep=",")
        if line[0] == chromosome:
            if line[21] == chromosome[3:] + band:
                if int(line[20]) in node_id_list:
                    populations_list[int(line[20])] = line[1]
                    tree_index_list[int(line[20])] = int(line[19])
                    chunk_index_list[int(line[20])] = int(line[18])

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

def get_bitset(n_set, samp_to_ind):
    e = [0] * int(len(samp_to_ind) / 2)
    for n in n_set:
        e[samp_to_ind[n]] += 1
    return e

def get_chunk(pos, ch):
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

def draw_tree(pos, p_list, n_list, t_list, c_list, ch, y_axis=False):

    print("Getting metadata...")
    poplabels_loc = "trees/relate_1000GP/data"
    md, populations, groups = argbd.simulations.read_in_pop_metadata(
        poplabels_loc + "/1000GP_Phase3.poplabels"
    )
    print("Groups and populations:")
    print(groups, populations)
    samp_to_ind_all = {}
    inds = {}
    i = 0
    for key, val in md.items():
        if val["ID"] not in inds:
            inds[val["ID"]] = i
            i += 1
        samp_to_ind_all[key] = inds[val["ID"]]
    samples_list = []
    bitset_list = []

    for node_id in n_list:
        population = p_list[node_id]
        tree_index = t_list[node_id]
        chunk_index = c_list[node_id]
        mapp = pickle.load(
            open(
                "trees/relate_1000GP/trees/chunks/"
                + ch
                + "/1000GP_Phase3_mask_prene_"
                + ch
                + "_"
                + population
                + "_map.pickle",
                "rb",
            )
        )
        ts = tszip.decompress(
            "trees/relate_1000GP/trees/chunks/"
            + ch
            + "/1000GP_Phase3_mask_prene_"
            + ch
            + "_chunk"
            + str(chunk_index)
            + "_"
            + population
            + ".trees.tsz"
        )
        t = ts.at_index(tree_index)
        ss = {np.where(mapp == s)[0][0] for s in t.samples(node_id)}
        bitset_all = [0] * int(len(samp_to_ind_all) / 2)
        for s in ss:
            bitset_all[samp_to_ind_all[s]] += 1
        bitset_list.append(bitset_all)
        samples_list.append(ss)

    styles = []
    chunk = get_chunk(pos, chromosome)
    ts = tszip.decompress(
        "trees/relate_1000GP/trees/chunks/"
        + ch
        + "/1000GP_Phase3_mask_prene_"
        + ch
        + "_chunk"
        + str(chunk)
        + ".trees.tsz"
    )
    t = ts.at(pos)
    print("chunk", chunk)
    scale = 10
    recorded = []
    for samples in samples_list:
        for ss in samples:
            styles.append(
                ".n"
                + str(ss)
                + " > .sym {transform: scale(0.25, "
                + str(scale)
                + ") translate(0, 0.8px); fill: "
                + colours[md[ss]["group"]]
                + "}"
            )
            recorded.append(ss)
    for ss in t.samples():
        if ss not in recorded:
            styles.append(
                ".n"
                + str(ss)
                + " > .sym {transform: scale(0.25, "
                + str(scale)
                + ") translate(0, 2px); fill: "
                + colours[md[ss]["group"]]
                + "}"
            )
        if ind_to_highlight != "None":
            if md[ss]["ID"] == ind_to_highlight:
                print("Highlighting sample", ss)
                print(md[ss])
                styles.append(
                    ".n"
                    +str(ss)
                    +" > .sym {transform: scale(0.25, "
                    +str(3*scale)
                    +") translate(0, 0.1px); fill: black"
                    +"}"
                )
    print(t.index)

    if node_to_calculate != "None":
        samples = [s for s in t.samples(node_to_calculate)]
        counts = {"EUR": 0, "AFR": 0, "SAS": 0, "EAS": 0, "AMR": 0}
        for s in samples:
            counts[md[s]["group"]] += 1
        counts_out = chromosome + " " + band + " " + str(pos) + " " + str(node_to_calculate) + " "
        for k, v in counts.items():
            print(k, v)
            counts_out += str(v/num_samples[k]) + " "
        print(counts_out)
    else:
        y_ticks = None
        if y_axis:
            y_ticks = [i * 10000 for i in range(0, 100)]
        svg = t.draw_svg(
            path=name + str(int(pos)) + "_.svg",
            # size=(250, 1000),
            size=(250,400),
            # size=(3000,1000),
            style="".join(styles),
            # node_labels={n:str(n) for n in t.nodes() if t.time(n) >= 5000},
            node_labels = {},
            symbol_size=1,
            mutation_labels={},
            # max_time="ts",
            y_axis=y_axis,
            y_ticks=y_ticks,
            y_label="    ",
        )
        # svg2pdf(
        #     bytestring=svg,
        #     write_to=name + str(int(pos)) + "_.pdf",
        #     dpi=600,
        # )

for p in positions:
    draw_tree(p, populations_list, node_id_list, tree_index_list, chunk_index_list, chromosome, y_axis=True)
