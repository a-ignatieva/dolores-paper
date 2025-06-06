import sys
import os
import glob

sys.path.append("..")

results_dir = "files/"
output_dir = "t2t/plots/"
index_file = results_dir + "children-to-parents.csv"
summary_file = results_dir + "Supplementary_Table_S1_significant_clades_details.csv"
output_file = output_dir + "counts.csv"

directory = os.fsencode(results_dir)
filelist = glob.glob(os.path.join(results_dir, '*carriers_.csv'))

par_to_ch = {}
with open(index_file, "r") as ifile:
    for line in ifile:
        line = line.strip().split(",")
        par_to_ch[line[3]] = (line[1], line[2])
print(par_to_ch)

frequencies = {}
with open(summary_file, "r") as ifile:
    for l, line in enumerate(ifile):
        if l > 1:
            line = line.strip().split(",")
            print(line[2])
            frequencies[line[2]] = (float(line[31]), float(line[32]), float(line[33]), float(line[34]), float(line[35]))

num_samples = {
    "EUR": 1006,
    "AFR": 1322,
    "SAS": 978,
    "EAS": 1008,
    "AMR": 694,
}

done = set()
redo = set()

with open(output_file, "w") as ofile:
    for filename in sorted(filelist):
        ofile.write(" " + ",")
        with open(filename, "r") as file:
            for line in file:
                line = line.strip().split(",")
                if line[0] in par_to_ch:
                    ofile.write(line[0] + ",")
        break
    ofile.write("\n")
with open(output_file, "a") as ofile:
    for filename in sorted(filelist):
        ofile.write(" " + ",")
        with open(filename, "r") as file:
            for line in file:
                line = line.strip().split(",")
                if line[0] in par_to_ch:
                    ofile.write(par_to_ch[line[0]][0] + "." + par_to_ch[line[0]][1] + ",")
        break
    ofile.write("\n")
with open(output_file, "a") as ofile:
    for filename in sorted(filelist):
        ofile.write(" " + ",")
        with open(filename, "r") as file:
            for line in file:
                line = line.strip().split(",")
                if line[0] in par_to_ch:
                    ofile.write(line[1] + ",")
        break
    ofile.write("\n")

with open(output_file, "a") as ofile:
    for filename in sorted(filelist):
        sv = filename.split("/")[-1]
        sv = sv.split("_")[0][3:]
        print(sv)
        done.add(sv)
        pop_counts = {"EUR": 0, "AFR": 0, "SAS": 0, "EAS": 0, "AMR": 0}
        hets = homs0 = homs1 = 0
        ofile.write(sv + ",")
        with open(filename, "r") as file:
            for line in file:
                line = line.strip().split(",")
                pop_counts[line[1]] += int(line[3])
                if line[0] in par_to_ch:
                    if int(line[3]) == 0:
                        # print("hom0", line[0], par_to_ch[line[0]])
                        homs0 += 1
                        ofile.write(str(0) + ",")
                    if int(line[3]) == 1:
                        # print("het", line[0], par_to_ch[line[0]])
                        hets += 1
                        ofile.write(str(0.5) + ",")
                    elif int(line[3]) == 2:
                        # print("hom1", line[0], par_to_ch[line[0]])
                        homs1 += 1
                        ofile.write(str(1) + ",")
            print(homs0, hets, homs1)
            for i, (k, v) in enumerate(pop_counts.items()):
                print(k, round(v/num_samples[k],2), frequencies[sv][i], abs(frequencies[sv][i] - round(v/num_samples[k],2)))
        ofile.write("\n")
        print("="*50)

for k in frequencies.keys():
    if k not in done:
        print(k)




