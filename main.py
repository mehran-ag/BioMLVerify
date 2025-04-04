from biomodel import *
import os


os.system("cls" if os.name == "nt" else "clear")

biomodel = BioModel()

biomodel.read_file("/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000315.xml")

print("\n\n")

# print(biomodel.getStoichiometricMatrix())

# print("\n\n")

# print(biomodel.getStoichiometricColumnNamesIndices())

# print("\n\n")

# print(biomodel.getStoichiometricRowNamesIndices())

# print("\n\n")

# biomodel.getElementInformationInStoichiometricMatrix(0,0)
biomodel.getThermoConversionMatrix()