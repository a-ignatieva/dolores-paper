#!/usr/bin/env python
# coding: utf-8

import stdpopsim
import tskit
import sys
import os

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

# python -m tgp_simulations_scr chr20 ~/Data/mask 1 25 5

chromosome = sys.argv[1]  # e.g. chr21
mask = sys.argv[2]  # location of genomic masks to use with Relate
rs = sys.argv[3]  # random seed
memory = int(sys.argv[4])  # total memory (will be divided by threads)
threads = int(sys.argv[5])  # number of threads for Relate

if mask == "None":
    mask = None
else:
    mask = mask + "/20140520." + chromosome + ".pilot_mask.fasta.gz"

relate_loc = "software/relate-master"
relate_lib_loc = "software/relate_lib-master"
simulation_loc = "trees/relate_1000GP_sim/trees/"
rec_map_file = "genetic_maps/HapMapII_GRCh37/genetic_map_GRCh37_" + chromosome + ".txt"

species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("AmericanAdmixture_4B11")
contig = species.get_contig(chromosome=chromosome, genetic_map="HapMapII_GRCh37", mutation_rate=model.mutation_rate)
mutation_rate = contig.mutation_rate
recombination_map = contig.recombination_map
samples = {"AFR": int(1322/2), "EUR": int(1006/2), "ASIA": int(1008/2), "ADMIX": int(694/2)}
Ne = species.population_size  # diploids

if os.path.exists(simulation_loc + str(rs) + "_sim_" + chromosome + ".trees"):
    print("Already have simulated trees, running Relate only")
    ts = tskit.load(simulation_loc + str(rs) + "sim_" + chromosome + ".trees")
    if ts.sequence_length > chromosome_lengths[chromosome]:
        ts = ts.keep_intervals([(0, chromosome_lengths[chromosome])])
        ts.dump(simulation_loc + str(rs) + "_sim_" + chromosome + ".trees")
else:
    print("Simulating " + chromosome, flush=True)
    print("mutation rate " +str(mutation_rate), flush=True)
    print("recombination rate " + str(recombination_map.mean_rate), flush=True)
    print("samples " + str(samples), flush=True)
    ts, rs = dolores.simulations.simulate_data(
        contig,
        model,
        samples,
    )
    print("random seed " + str(rs), flush=True)
    print(
        ts.num_trees,
        ts.num_samples,
        ts.first().interval,
        ts.last().interval,
        flush=True,
    )
    ts.dump(
        simulation_loc + str(rs) + "_sim_" + chromosome + ".trees"
    )

    # dolores.simulations.write_poplabels(ts, simulation_loc + "simulated)


print("Running Relate", flush=True)
dolores.simulations.run_relate(
    ts,
    rec_map_file,
    mutation_rate,
    2 * Ne,
    relate_loc,
    mask=mask,
    memory=int(memory / threads),
    consistency=False,
    postprocess=False,
    randomise=False,
    threads=threads,
    quiet=False,
)
os.system(
    "mv relate.anc "
    + simulation_loc
    + "relate_sim_"
    + chromosome
    + "_part0.anc"
)
os.system(
    "mv relate.mut "
    + simulation_loc
    + "relate_sim_"
    + chromosome
    + "_part0.mut"
)
os.system(
    "gzip "
    + simulation_loc
    + "relate_sim_"
    + chromosome
    + "_part0.anc"
)
os.system(
    "gzip "
    + simulation_loc
    + "relate_sim_"
    + chromosome
    + "_part0.mut"
)
os.system(
    relate_loc
    + "/bin/RelateFileFormats --mode ConvertToTreeSequence -i " + simulation_loc + "relate_sim_"
    + chromosome
    + "_part0 -o " + simulation_loc + "relate_sim_"
    + chromosome
    + "_part0;"
)

print("Running final cleanup", flush=True)
dolores.simulations.clean_relate()
print("Done", flush=True)
