import sys
import os
import gzip
import glob
import numpy as np
import scipy as sp
import tqdm

sys.path.append("..")

extract_kmers = False

results_dir = "files/"
data_dir = "t2t/genomes/"
ann_dir = "t2t/annotations/"
output_dir = "t2t/kmers/"
band_file = "files/bands_to_flanking_genes_.csv"
directory = os.fsencode(output_dir)
k = 20

dna_dict = {"A": "T", "C": "G", "T": "A", "G": "C"}

def reverse_seq(string):
    new_string = ""
    for s in string:
        new_string += dna_dict[s]
    return new_string

def flip_seq(string):
    new_string = string[::-1]
    return new_string

with open("files/counts_HPRC_updated.csv", "r") as file:
    j = 0
    IDs = {}
    carriers = {}
    for l, line in enumerate(file):
        line = line.strip().split(",")
        if l == 0:
            for ll in line:
                if ll != "":
                    IDs[ll] = j
                    j += 1
        elif l > 1:
            carriers[line[0]] = [int(ll) for ll in line[1:] if ll != ""]

with open(band_file, "r") as band_file:
    for lline in band_file:
        lline = lline.strip().split(",")
        band = lline[0]
        startgene = lline[1]
        endgene = lline[2]

        startgene_dir = endgene_dir = None
        output_file = output_dir + band + "_kmer_counts_"
        summary_file = output_dir + band + "_kmer_summary.csv.gz"
        genes = {startgene, endgene}

        if extract_kmers:
            directory = os.fsencode(results_dir)
            filelist = glob.glob(os.path.join(data_dir, '*fa.gz'))
            # filename = "Homo_sapiens-GCA_018466835.1-softmasked.fa.gz"
            # filename = "Homo_sapiens-GCA_018469695.1-softmasked.fa.gz"
            # filename = "Homo_sapiens-GCA_018471545.1-softmasked.fa.gz"

            for f, filename in enumerate(sorted(filelist)):
                print(filename)

                genefile = filename.split("/")[-1]
                genefile = genefile.split("-")
                genefile = genefile[0] + "-" + genefile[1]
                if os.path.exists(ann_dir + genefile + "-2022_07-genes.gtf.gz"):
                    genefile = genefile + "-2022_07-genes.gtf.gz"
                else:
                    genefile = genefile + "-2022_08-genes.gtf.gz"

                coords = {}
                with gzip.open(ann_dir + genefile, "rt") as gfile:
                    print("Reading genes", genefile)
                    for l, line in enumerate(gfile):
                        line = line.strip()
                        if l == 0:
                            line = line.split(" ")
                            ref = line[1]
                        if l > 4:
                            line = line.split("\t")
                            if line[2] == "gene":
                                r = line[0]
                                s = int(line[3])
                                e = int(line[4])
                                st = line[6]
                                details_ = line[-1].split(";")
                                details = {}
                                for d in details_:
                                    d = [dd.strip() for dd in d.split('"')]
                                    if len(d) > 1 and d[0] not in details:
                                        details[d[0]] = d[1]
                                if "gene_name" in details:
                                    gene_name = details["gene_name"]
                                    if gene_name in genes:
                                        print(gene_name, r)
                                        if gene_name in coords:
                                            print("Gene on two contigs", coords[gene_name][0], r)
                                        else:
                                            coords[gene_name] = [r, s, e, st]
                                        if f == 0:
                                            if gene_name == startgene:
                                                startgene_dir = st
                                            else:
                                                endgene_dir = st
                if coords[startgene][0] == coords[endgene][0]:
                    contig = [coords[startgene][0]]
                    start = [min(coords[startgene][1], coords[endgene][1])]
                    end = [max(coords[startgene][2], coords[endgene][2])]
                    strand = [coords[startgene][3]]
                    rev = [coords[startgene][3] != startgene_dir]
                else:
                    print("---> Different contigs")
                    contig = [coords[startgene][0], coords[endgene][0]]
                    rev = [coords[startgene][3] != startgene_dir, coords[endgene][3] != endgene_dir]
                    strand = [coords[startgene][3], coords[endgene][3]]
                    start = [None, None]
                    end = [None, None]
                    for i, (gene, gene_dir) in enumerate([(startgene, startgene_dir), (endgene, endgene_dir)]):
                        if (i == 0 and coords[gene][3] == gene_dir) or (i == 1 and coords[gene][3] != gene_dir):
                            start[i] = coords[gene][1]
                            end[i] = 1000000000
                        else:
                            start[i] = 0
                            end[i] = coords[gene][2]
                print(contig, start, end, strand, rev)

                chunk = ""
                for i in range(len(contig)):
                    sequence = ""
                    readnext = False
                    with gzip.open(filename, "rt") as file:
                        for line in file:
                            line = line.strip()
                            if line[0] != ">":
                                if readnext:
                                    sequence += line
                            elif not readnext:
                                line = line[1:].split()
                                name = line[0]
                                if name == contig[i]:
                                    line_ = line[2].split(".")
                                    if line_[0][0] == "p":
                                        line_ = line_[0].split(":")
                                        ind = line_[1]
                                        mp = "none"
                                    else:
                                        ind = line_[0]
                                        mp = line_[2]
                                    print(name, ind, mp)
                                    readnext = True
                            else:
                                break
                    print("Contig length:", len(sequence))
                    if rev[i]:
                        chunk_ = sequence[max(0, start[i] - 2) : min(end[i], len(sequence))]
                        chunk_ = flip_seq(chunk_)
                        chunk_ = reverse_seq(chunk_)
                        print("Reversing...")
                    else:
                        chunk_ = sequence[max(0, start[i] - 1) : min(end[i] + 1, len(sequence))]
                    chunk += chunk_

                print(chunk[0:20] + "..." + chunk[-21:-1])
                print("Search length:", len(chunk))
                print("Total kmers:", len(chunk) - k + 1)

                kmer_counts = {}
                kmer_pos = {}
                k_check = 0
                for i in range(len(chunk) - k + 1):
                    kmer = chunk[i:(i+k)]
                    if kmer in kmer_counts:
                        kmer_counts[kmer] += 1
                        kmer_pos[kmer].append(i)
                    else:
                        kmer_counts[kmer] = 1
                        kmer_pos[kmer] = [i]
                    k_check += 1
                # for k, v in kmer_counts.items():
                #     print(k + " : " + str(v))
                # print("Max kmer count:", max([v for v in kmer_counts.values()]))
                # kmer_dist = {i: 0 for i in range(1 + max([v for v in kmer_counts.values()]))}
                # for k, v in kmer_counts.items():
                #     kmer_dist[v] += 1
                # for k, v in kmer_dist.items():
                #     print(str(k) + " : " + str(v))
                print("Total kmers:",  k_check)

                print("Writing to file...")
                with gzip.open(output_file + ind + "." + mp + ".csv.gz", "wt") as ofile:
                    for kmer in kmer_counts.keys():
                        ofile.write(kmer + "," + str(kmer_counts[kmer]) + "," + ";".join(str(p) for p in kmer_pos[kmer]) + "\n")
                print("="*50)


        filelist = glob.glob(os.path.join(output_dir, band + '_kmer_counts_*.csv.gz'))
        kmers_ref = set()
        kmers_ref_list = {}

        with gzip.open(output_dir + band + "_kmer_counts_T2T-CHM13v2.none.csv.gz", "rt") as file:
            for line in file:
                line = line.strip().split(",")
                kmers_ref.add(line[0])
                # if int(line[1]) == 1:
                #     kmers_ref.add(line[0])
        total_kmers = len(kmers_ref)
        # print("Total kmers:", total_kmers)
        with gzip.open(output_dir + band + "_kmer_counts_T2T-CHM13v2.none.csv.gz", "rt") as file:
            for line in file:
                line = line.strip().split(",")
                if int(line[1]) == 1:
                    kmer = line[0]
                    kmer_ = flip_seq(reverse_seq(kmer))
                    if kmer_ not in kmers_ref:
                        kmers_ref_list[kmer] = int(line[2])
        kmers_ref = set(kmer for kmer in kmers_ref_list)
        unique_kmers = len(kmers_ref)
        # print("Total unique kmers:", unique_kmers)

        inv_count = np.zeros(len(filelist))
        dup_count = np.zeros(len(filelist))
        del_count = np.zeros(len(filelist))
        normal_count = np.zeros(len(filelist))
        status = np.zeros(len(filelist))
        with gzip.open(summary_file, "wt") as sfile:
            for i, filename in enumerate(sorted(filelist)):
                kmers_seq = set()
                filename_ = filename.strip().split("_")[3].split(".")
                name = filename_[0]
                mp = filename_[1]
                with gzip.open(filename, "rt") as file:
                    for line in file:
                        line = line.strip().split(",")
                        kmer = line[0]
                        kmers_seq.add(kmer)
                        ct = int(line[1])
                        if ct == 1:
                            if kmer not in kmers_ref:
                                kmer_ = flip_seq(reverse_seq(kmer))
                                if kmer_ in kmers_ref:
                                    inv_count[i] += 1
                            else:
                                normal_count[i] += 1
                        if kmer in kmers_ref and ct > 1:
                            dup_count[i] += 1
                    del_count[i] = sum(1 for kmer in kmers_ref if kmer not in kmers_seq and flip_seq(reverse_seq(kmer)) not in kmers_seq)

                out = band + "," + str(total_kmers) + "," + str(unique_kmers) + "," + name + "." + mp + ","
                out += str(len(kmers_seq)) + ","
                out += str(normal_count[i]) + "," + str(inv_count[i]) + "," + str(dup_count[i]) + "," + str(del_count[i]) + ","
                if name + "." + mp in IDs:
                    out += str(carriers[band][IDs[name + "." + mp]])
                    status[i] = carriers[band][IDs[name + "." + mp]]
                else:
                    out += "NA"
                    status[i] = -1
                # print(out)
                out += "\n"
                sfile.write(out)

        inv_count = [inv_count[i] for i in range(len(status)) if status[i] != -1]
        dup_count = [dup_count[i] for i in range(len(status)) if status[i] != -1]
        del_count = [del_count[i] for i in range(len(status)) if status[i] != -1]
        normal_count = [normal_count[i] for i in range(len(status)) if status[i] != -1]
        status = [status[i] for i in range(len(status)) if status[i] != -1]

        out = band
        for sv, A in [("inversion", inv_count), ("duplication", dup_count), ("deletion", del_count)]:
            # print(sv + " counts:")
            for p in [0, 1]:
                pp = np.mean([A[i] for i in range(len(status)) if status[i] == p])
                # print(p, pp)
                out += " " + str(pp)
            pp = sp.stats.pearsonr(A, status)
            # print(pp[0])
            out += " " + str(pp[0])
        print(out)




