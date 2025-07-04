from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path, PurePath



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


def verify_model(file_path, file_name):

    # Ensure file_path is a string or Path
    if isinstance(file_path, (Path, PurePath)):
        file_path = str(file_path)
    elif not isinstance(file_path, str):
        raise TypeError('The file path should be a string or a Path object.')

    full_path = os.path.join(file_path, file_name)

    if not os.path.isfile(full_path):
        raise FileNotFoundError(f'Model source file `{full_path}` does not exist.')
    
    biomodel = BioModel()

    biomodel.read_file(full_path)

    if not biomodel.checkMassActionKinetics():

        utility.message_printer(f"\nModel {file_name} has (a) reaction(s) not governed by \"Mass Action\" kinetics and\nis not eligible for verification check\n\n\n", color='red', style='normal')

        return
    
    if not biomodel.checkModelReversibility():

        utility.message_printer(f"\nModel {file_name} has (a) irrversible reaction(s) and\nis not eligible for verification check\n\n\n", color='red', style='normal')

        return
    
    biomodel.KineticConstantsThermoCompatibilty("on")



