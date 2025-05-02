from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from functions import *

os.system("cls" if os.name == "nt" else "clear")

path = "/Users/makb047/UoA/Codes/Mass_Actions"

# read_files_in_folder(path)


biomodel = BioModel()

biomodel.read_file("/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000002.xml")

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