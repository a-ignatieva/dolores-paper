#!/usr/bin/env python
# coding: utf-8

import sys
import os
import random
import tskit
import msprime
import pyslim
import numpy as np

sys.path.append("..")
import dolores.simulations

def write_poplabels(ts, filename):
    """
    Write a file with the population labels for each sample
    :param ts: simulated ts
    :param filename: where to write .poplabels file (no extension)
    :return:
    """
    with open(filename + ".poplabels", "w") as file:
        file.write("sample population group sex\n")
        for s in ts.individuals():
            node = ts.node(s.nodes[0])
            pop = ts.population(node.population).metadata["name"]
            file.write("tsk_" + str(s.id) + " " + pop + " " + pop + " NA\n")
            if len(s.nodes) > 1:
                node = ts.node(s.nodes[1])
                pop_ = ts.population(node.population).metadata["name"]
                if pop_ != pop:
                    sys.exit(
                        "Warning: can't have different population IDs for the same individual."
                    )

invlength = int(sys.argv[1])
samplesize = int(sys.argv[2])  # diploids
run = int(sys.argv[3])

outname =  "slim_" + str(invlength) + "_" + str(samplesize) + "_" + str(run)

relate_loc = "/Users/ignatiev/Desktop/relate"

# Simulate trees using SLiM script
# These will be called slim.trees
stop = False
while not stop:
    if invlength > 0:
        print("Simulating with inversion")
        print("slim -d invlength=" + str(invlength) + " slim-sim-with-inversion.slim")
        os.system("slim -d invlength=" + str(invlength) + " slim-sim-with-inversion.slim")
    else:
        print("Simulating without inversion")
        os.system("slim slim-sim-neutral.slim")

    ts = tskit.load("slim.trees")
    print(ts.num_trees)
    print(ts.num_mutations)
    if invlength > 0 and ts.num_mutations == 2:
        stop = True

Ne = 10000  # diploids
mutation_rate = 1e-8
recombination_rate = 1e-8

demography = msprime.Demography.from_tree_sequence(ts)
for pop in demography.populations:
    if pop.name == 'p1':
        pop.initial_size = Ne
recap = pyslim.recapitate(
    ts, 
    demography=demography,
    recombination_rate=recombination_rate, 
    random_seed=1,
)

seed = random.randrange(100001)
print("Random seed:", seed)
rng = np.random.default_rng(seed = seed)
alive_inds = pyslim.individuals_alive_at(recap, 0)
keep_indivs = rng.choice(alive_inds, samplesize, replace=False)
keep_nodes = []
for i in keep_indivs:
    keep_nodes.extend(recap.individual(i).nodes)

recap = recap.simplify(keep_nodes, keep_input_roots=False)

samp_to_ind = {}
for s in recap.individuals():
    node = recap.node(s.nodes[0])
    for n in s.nodes:
        samp_to_ind[n] = s.id

if invlength > 0 and recap.num_mutations == 0:
    sys.exit("No inversions here!")
elif invlength > 0:
    print("Inversion is above nodes:")
    for m in recap.mutations():
        print(m.node)
        p = recap.site(m.site).position
        print("Inversion position:", p)
        t = recap.at(p)
        samples = [s for s in t.samples(m.node)]
        f = len(samples)/t.num_samples()
        print("Inversion frequency:", f)
        if 0.05 < f < 0.95:
            carriers = np.zeros(samplesize)
            for s in samples:
                carriers[samp_to_ind[s]] += 0.5
            with open("carriers.txt", "w") as file:
                for s in recap.individuals():
                    file.write("tsk_" + str(s.id) + "," + str(carriers[s.id]) + "\n")
    print("="*10)
else:
    print("No inversion simulated.")
    with open("carriers.txt", "w") as file:
        for s in recap.individuals():
            file.write("tsk_" + str(s.id) + ",0\n")

# WARNING: need keep=False (otherwise will have an issue with running Relate later)
recap = msprime.sim_mutations(recap, rate=mutation_rate, random_seed=11, discrete_genome=True, keep=False)
sites_to_delete = [s.id for s in recap.sites() if len(s.mutations) > 1]
print("Deleting homoplasic sites:", len(sites_to_delete), flush=True)
recap = recap.delete_sites(sites_to_delete)
recap.dump(outname + "_recap.trees")
print("Number of mutations:", recap.num_mutations)
write_poplabels(recap, outname + "_recap")

with open(outname + "_recap.popsize", "w") as file:
    file.write("0,20000\n10000000,NA\n")

recombination_map = msprime.RateMap(
    position=[0, ts.sequence_length],
    rate=[recombination_rate],
)

dolores.simulations.flat_recombination_map(recombination_rate, recap.sequence_length)

dolores.simulations.run_relate(
    recap,
    rec_map_file="dummy_map.txt",
    mutation_rate=mutation_rate,
    Ne=Ne*2,
    relate_loc=relate_loc,
    quiet=False,
)

os.system(relate_loc + '/bin/RelateFileFormats --mode ConvertToTreeSequence -i relate -o ' + outname)
os.system("mv " + outname + ".trees " + outname + "_relate.trees")
os.system("mv " + "slim.trees " + outname + ".trees")
os.system("cp " + outname + "_recap.poplabels " + outname + "_relate.poplabels")
os.system("cp " + outname + "_recap.popsize " + outname + "_relate.popsize")
