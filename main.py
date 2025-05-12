from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from functions import *



os.system("cls" if os.name == "nt" else "clear")

folder_path = "/Users/makb047/UoA/Codes/Mass_Actions"

# read_files_in_folder(folder_path)


# file_path = "/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000500.xml"

file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/huang_ferrell_1996_original.cellml"

biomodel = BioModel()

biomodel.read_file(file_path)

print("\n")

biomodel.getStoichiometricMatrix("on")

biomodel.getForwardStoichiometricMatrix("on")

biomodel.getReverseStoichiometricMatrix("on")

large_array = biomodel.getThermoConversionMatrix("on")

df = pd.DataFrame(large_array)

print("\nThe Kinetic to Thermodynamic Conversion Matrix is:\n",df)

biomodel.checkModelReversibility(printing="on")

biomodel.getKineticRateConstantsVector("on")

biomodel.KineticConstantsThermoCompatibilty("on")