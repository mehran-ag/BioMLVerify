from bioml import *
import os
import atexit


from _modules._utility import time_counter

from time import perf_counter as t; t0 = t()

os.system("cls" if os.name == "nt" else "clear")

atexit.register(lambda: print("\n" * 2))


folder_path = "SBML and CellML test models"

file_name = "BIOMD0000000085.xml"

bioml = BioML()

bioml.read_file(folder_path, file_name)

bioml.verify_model(mass_balance=False, charge_balance=False, printing=True)


time_counter(t, t0)  # This function calculates the execution time and prints on the screen