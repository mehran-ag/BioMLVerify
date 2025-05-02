from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



def read_files_in_folder(path):


    reversibles = []

    for filename in os.listdir(path):

        if not filename.endswith('.xml'): continue

        fullname = os.path.join(path, filename)

        biomodel = BioModel()

        biomodel.read_file(fullname)

        if biomodel.checkModelReversibility():
            print(f"\nModel {biomodel.file_name} is ALL REVERSIBLE")
            reversibles.append(biomodel.file_name)

    print(len(reversibles))

    for m in reversibles:

        print(f"\n{m}")


    with open("/Users/makb047/UoA/Codes/Mass_Actions/reversibles.txt", "w") as file:
        for item in reversibles:
            file.write(item + "\n")

    return 