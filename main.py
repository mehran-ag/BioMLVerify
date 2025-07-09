from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from functions import *



os.system("cls" if os.name == "nt" else "clear")

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs"

# file_name = "modified_huang_ferrell_1996.cellml"




file_path = "/Users/makb047/UoA/Codes/Mass_Actions"

file_name = "BIOMD0000000019.xml"

# verify_model(file_path, file_name)

verify_bunch_SBML_models(file_path)





# file_path = "/Users/makb047/UoA/Codes/aguda_b_1999/aguda_b_1999-original.cellml"

# file_path = "/Users/makb047/UoA/Codes/NitrosylBromide_BioML/NitrosylBromide-BioML.cellml"

# file_path = "/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000192.xml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/modified_huang_ferrell_1996.cellml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/reactions_set.cellml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/aguda_b_1999.cellml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/NitrosylBromide.cellml"

# biomodel = BioModel()

# biomodel.read_file(file_path)

# print("\n")

# biomodel.checkMassActionKinetics("on")

# biomodel.getStoichiometricMatrix("on")

# biomodel.getStoichiometricColumnNamesIndices("on")

# biomodel.getStoichiometricRowNamesIndices("on")

# biomodel.getForwardStoichiometricMatrix("on")

# biomodel.getReverseStoichiometricMatrix("on")

# # large_array = biomodel.getThermoConversionMatrix("on")

# # # df = pd.DataFrame(large_array)

# # # print("\nThe Kinetic to Thermodynamic Conversion Matrix is:\n",df)

# biomodel.checkModelReversibility(printing="on")

# biomodel.getKineticRateConstantsVector("on")

# biomodel.KineticConstantsThermoCompatibilty("on")

