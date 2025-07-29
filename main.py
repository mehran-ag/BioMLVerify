from bioml import *
import os
import atexit

from functions import *

from utility import printer

from time import perf_counter as t; t0 = t()



os.system("cls" if os.name == "nt" else "clear")

atexit.register(lambda: print("\n" * 2))





# folder_path = "/Users/makb047/UoA/Codes/NitrosylBromide_BioML"

# folder_path = "/Users/makb047/UoA/Codes/ATP Hydrolysis"

# folder_path = "/Users/makb047/UoA/Codes/SBML_Models"

file_name = "BIOMD0000000150.xml"

# file_name = "ATP_Hydrolysis.cellml"

# file_name = "modified_huang_ferrell_1996.cellml"

# file_name = "aguda_b_1999.cellml"

# file_name = "aguda_b_1999-original.cellml"

# file_name = "NitrosylBromide-BioML.cellml"


# verify_bunch_models(folder_path)

folder_path = "/Users/makb047/UoA/Codes/Mass_Actions"


bioml = BioML()

bioml.read_file(folder_path, file_name)

bioml.verify_model(mass_balance=False, charge_balance=False, printing=True)





elapsed = t() - t0

elapsed_str = (
    f"{int(elapsed // 60)} minute{'s' if int(elapsed // 60) != 1 else ''} and "
    f"{int(elapsed % 60)} second{'s' if int(elapsed % 60) != 1 else ''}"
) if elapsed >= 60 else f"{int(elapsed)} second{'s' if int(elapsed) != 1 else ''}"

printer("\n\nExecution time: ", elapsed_str, text_color='light_yellow', text_style='dim')