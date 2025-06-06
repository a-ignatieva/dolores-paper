import os
import numpy as np
from sklearn.metrics import adjusted_rand_score

loc = "/Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/"


for il, invlength in enumerate([50000, 100000, 200000]):
    for i in range(100):
        print("="*80)
        print(i)
        print("=" * 80)
        nwindows = int(5000000/10000)
        loc_ = loc + "slim_" + str(invlength) + "_50_" + str(i) + "/"
        if os.path.exists(loc_ + "dolores_breakpoints.txt") and os.path.exists(loc_ + "dolores_prediction.txt"):
            filename = "/Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/slim_" + str(invlength) + "_50_" + str(i) + "/"
            filename_ = filename + "slim_" + str(invlength) + "_50_" + str(i) + "_asaph"
            os.system("cd " + loc_)
            os.system("asaph_pca --workdir " + filename_ + " pca --vcf " + filename + "samples.vcf --min-inversion-fraction 0.01")
            os.system("asaph_localize --workdir " + filename_ + " association-tests --components 1 2 3 4 --vcf " + filename + "samples.vcf")
            os.system("asaph_localize --workdir " + filename_ + " detect-boundaries --component 1 --n-windows " + str(
                nwindows) + " > " + filename + "asaph_breakpoints.txt")
            os.system("asaph_localize --workdir " + filename_ + " detect-boundaries --component 2 --n-windows " + str(
                nwindows) + " >> " + filename + "asaph_breakpoints.txt")
            os.system("asaph_localize --workdir " + filename_ + " detect-boundaries --component 3 --n-windows " + str(
                nwindows) + " >> " + filename + "asaph_breakpoints.txt")
            os.system("asaph_localize --workdir " + filename_ + " detect-boundaries --component 4 --n-windows " + str(
                nwindows) + " >> " + filename + "asaph_breakpoints.txt")
            os.system("asaph_genotype cluster --workdir " + filename_ + " --components 1 --n-clusters 3 --predicted-labels-fl " + filename + "asaph_prediction.txt")
            os.system("cat " + filename + "asaph_breakpoints.txt")

            carriers = np.zeros(50)
            with open(filename + "carriers.txt", "r") as file:
                for i, line in enumerate(file):
                    line = line.strip().split(",")
                    carriers[i] = int(float(line[1])*2)
            predictions = np.zeros(50)
            with open(filename + "asaph_prediction.txt", "r") as file:
                for line in file:
                    line = line.strip().split(",")
                    label = -1
                    for j, n in enumerate(line):
                        if j == 0:
                            label = int(n)
                        else:
                            n = int(n.split("_")[1])
                            predictions[n] = label
            print(predictions)
            print(carriers)
            print(adjusted_rand_score(carriers, predictions))
            with open(filename + "asaph_score.txt", "w") as file:
                file.write(str(adjusted_rand_score(carriers, predictions)) + "\n")

