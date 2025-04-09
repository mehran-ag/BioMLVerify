from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.system("cls" if os.name == "nt" else "clear")

biomodel = BioModel()

biomodel.read_file("/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000315.xml")

print("\n\n")

biomodel.getForwardStoichiometricMatrix(printing="on")

biomodel.getReverseStoichiometricMatrix(printing="on")

biomodel.getStoichiometricMatrix(printing="on")

biomodel.getElementInformationInStoichiometricMatrix(1,5)

large_array = biomodel.getThermoConversionMatrix()

df = pd.DataFrame(large_array)

print("\nThe Kinetic to Thermodynamic Conversion Matrix is:\n",df)

biomodel.getKineticRateConstantsVector("on")

check = biomodel.KineticConstantsThermoCompatibilty("on")