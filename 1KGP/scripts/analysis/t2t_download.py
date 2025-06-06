import os

with open("ID_list.txt", "r") as ifile:
    for ID in ifile:
        ID = ID.strip()
        ofile = "Homo_sapiens-" + ID + "-softmasked.fa.gz"
        if os.path.exists(ofile):
            print(ID + " already downloaded")
        else:
            print(ID + "downloading...")
            url = "https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/" + ID + "/ensembl/genome/Homo_sapiens-" + ID + "-softmasked.fa.gz"
            os.system("curl -o " + ofile + " " + url)



