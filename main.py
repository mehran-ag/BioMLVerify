from biomodel import *
import os
import atexit

from functions import *

from utility import printer

from time import perf_counter as t; t0 = t()



os.system("cls" if os.name == "nt" else "clear")

atexit.register(lambda: print("\n" * 2))

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs"

# file_name = "modified_huang_ferrell_1996.cellml"




folder_path = "/Users/makb047/UoA/Codes/SBML_Models"

file_name = "BIOMD0000000607.xml"

verify_model(folder_path, file_name)

# verify_bunch_SBML_models(folder_path)





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

printer("\n\nElapsed: ", f"{t() - t0:.4f} Seconds", text_color='light_yellow', text_style='dim')