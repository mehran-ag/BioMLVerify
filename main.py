from bioml import *
import os
import atexit

from functions import *

from utility import printer

from time import perf_counter as t; t0 = t()



os.system("cls" if os.name == "nt" else "clear")

atexit.register(lambda: print("\n" * 2))

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs"

# file_name = "modified_huang_ferrell_1996.cellml"




# folder_path = "/Users/makb047/UoA/Codes/SBML_Models"

# file_name = "BIOMD0000000283.xml"

# verify_model(folder_path, file_name)

# verify_bunch_SBML_models(folder_path)



# file_path = "/Users/makb047/UoA/Codes/aguda_b_1999/aguda_b_1999-original.cellml"

file_path = "/Users/makb047/UoA/Codes/NitrosylBromide_BioML/NitrosylBromide-BioML.cellml"

# file_path = "/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000192.xml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/modified_huang_ferrell_1996.cellml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/reactions_set.cellml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/aguda_b_1999.cellml"

# file_path = "/Users/makb047/UoA/Codes/CellML_Model_Verification/docs/NitrosylBromide.cellml"

bioml = BioML()

bioml.read_file(file_path)

# print("\n")

# bioml.check_mass_action_kinetics("on")

# bioml.get_stoichiometric_matrix("on")

# bioml.get_stoichiometric_column_names_indices("on")

# bioml.get_stoichiometric_row_names_indices("on")

# bioml.get_forward_stoichiometric_matrix("on")

# biomlm.get_reverse_stoichiometric_matrix("on")

# large_array = biomlmodel.get_thermo_conversion_matrix("on")

# df = pd.DataFrame(large_array)

# print("\nThe Kinetic to Thermodynamic Conversion Matrix is:\n",df)

# biomlm.check_model_reversibility(printing="on")

# bioml.get_kinetic_rate_constants_vector("on")

# bioml.check_kinetic_constants_thermo_compatibility("on")

elapsed = t() - t0

elapsed_str = (
    f"{int(elapsed // 60)} minute{'s' if int(elapsed // 60) != 1 else ''} and "
    f"{int(elapsed % 60)} second{'s' if int(elapsed % 60) != 1 else ''}"
) if elapsed >= 60 else f"{int(elapsed)} second{'s' if int(elapsed) != 1 else ''}"

printer("\n\nElapsed time: ", elapsed_str, text_color='light_yellow', text_style='dim')