from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path, PurePath

import gc



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

        utility.message_printer(f"\nModel {file_name} has (a) reaction(s) not governed by \"Mass Action\" kinetics and\nis not eligible for verification check\n\n\n", color='red')

        return
    
    if not biomodel.checkModelReversibility():

        utility.message_printer(f"\nModel {file_name} has (a) irrversible reaction(s) and\nis not eligible for verification check\n\n\n", color='red')

        return
    
    biomodel.KineticConstantsThermoCompatibility("on")
    del biomodel





def verify_bunch_SBML_models(folder_path):

    check_results = pd.DataFrame(columns=["Model Name", "Mass Action", "Reversible", "Plausible"])



    for file_name in os.listdir(folder_path):

        if not file_name.endswith('.xml'): continue

        full_path = os.path.join(folder_path, file_name)

        if not os.path.isfile(full_path):
            raise FileNotFoundError(f'Model source file `{full_path}` does not exist.')
        
        
        biomodel = BioModel()

        biomodel.read_file(full_path)

        mass_action: bool = None

        reversible: bool = None

        plausible: bool = None

        error: str = None

        try:


            if biomodel.checkMassActionKinetics():

                mass_action = True

                if biomodel.checkModelReversibility():

                    reversible = True

                    if biomodel.KineticConstantsThermoCompatibility():

                        plausible = True

                    else:

                        plausible = False

                else:

                    reversible = False

            else:

                mass_action = False

        except Exception as e:

            error = str(e)

        check_results.loc[len(check_results)] = [file_name, mass_action, reversible, plausible]

    excel_full_path = os.path.join(folder_path, "check_results.xlsx")
        
    check_results.to_excel(excel_full_path, index=False)