from biomodel import *
import os


os.system("cls" if os.name == "nt" else "clear")

biomodel = BioModel()

biomodel.read_file("/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000315.xml")

print("\n\n")

biomodel.getForwardStoichiometricMatrix(printing="on")

biomodel.getReverseStoichiometricMatrix(printing="on")

biomodel.getStoichiometricMatrix(printing="on")