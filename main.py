from biomodel import *
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.system("cls" if os.name == "nt" else "clear")

path = "/Users/makb047/UoA/Codes/Mass_Actions"

for filename in os.listdir(path):

    if not filename.endswith('.xml'): continue

    fullname = os.path.join(path, filename)

    biomodel = BioModel()

    biomodel.read_file(fullname)

    if biomodel.checkModelReversibility():
        print(f"\nModel {biomodel.file_name} is ALL REVERSIBLE")

# biomodel = BioModel()

# biomodel.read_file("/Users/makb047/UoA/Codes/Mass_Actions/BIOMD0000000267.xml")

# print("\n")

# large_array = biomodel.getThermoConversionMatrix()

# df = pd.DataFrame(large_array)

# print("\nThe Kinetic to Thermodynamic Conversion Matrix is:\n",df)

# biomodel.checkModelReversibility(printing="on")