#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import sys

results_file = (
    "files/clades_pvalues_all_events_varNe.csv"
)
new_file = (
    "files/clades_pvalues_all_events_annotated_varNe.csv"
)
bands_file = "files/cytoBand.txt"

bands = defaultdict(list)
with open(bands_file, "r") as file:
    for line in file:
        line = line.split()
        chromosome = line[0]
        start = int(line[1])
        end = int(line[2])
        name = line[3]
        bands[chromosome].append((start, end, name))
# print(bands['chr1'])
# print(bands["chr1"])

with open(results_file, "r") as file:
    with open(new_file, "w") as new_file:
        for l, line in enumerate(file):
            if l == 0:
                new_file.write(line[:-1] + ",band\n")
            else:
                line_ = line.split(sep=",")
                chromosome = line_[0]
                ch = chromosome[3:]
                start = int(line_[14])
                end = int(line_[15])
                found = False
                # print(chromosome, bands[chromosome])
                for band in bands[chromosome]:
                    # print(band)
                    if band[0] <= start and end <= band[1]:
                        new_file.write(line[:-1] + "," + ch + band[2] + "\n")
                        found = True
                        break
                if not found:
                    for band in bands[chromosome]:
                        if band[0] <= (start + end)/2 <= band[1]:
                            new_file.write(line[:-1] + "," + ch + band[2] + "\n")
                            found = True
                            break
                if not found:
                    print(chromosome, start, end)
                    sys.exit("Can't find band!")

