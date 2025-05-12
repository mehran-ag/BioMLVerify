import libsbml
import sys
import os
import traceback
import matrix_constructor
import model_checker
import exceptions
import utility
import warnings
import time
import sympy as sp

import sbml_reader

import cellml_reader

from classes.cReaction import *
from classes.cModel import *
from classes.cSpecies import *
from classes.cParameter import *
from classes.cSpeciesReference import *
from classes.cFunctionDefinition import *


class BioModel(object):
    '''
    This class creates a model, either CellML or SBML, that will be processed by other classes
    '''

    def __init__(self):
        '''
        Initializes the ModelReader by reading path for a file or a container of files
        '''
        
        self._file_path = None
        self._file_name= None
        self._file_format = None
        self._biomodel = None
        self._matrix_constructor = matrix_constructor.MatrixConstructor()
        self._model_checker = model_checker.ModelChecker()
        self._sbml_reader = sbml_reader.SbmlReader()
        self._cellml_reader = cellml_reader.CellmlReader()


    @property
    def file_name(self):
        return self._file_name


    def read_file(self, file_path):

        try:
            if not isinstance(file_path, str):
                raise TypeError("File path must be a string.")
            
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
            
            self._file_path = file_path
            self._file_name = os.path.basename(file_path)
            self._file_format = os.path.splitext(file_path)[1][1:]
            
        except TypeError as e:
            utility.message_printer("File not read!", color="red")
            utility.printer("\nERROR: ", e)
            return

        except FileNotFoundError as e:
            utility.message_printer("File not read!", color="red")
            utility.printer("\nERROR: ", e)
            return
        
        except Exception as e:
            utility.message_printer("File not read!", color="red")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return

        else:
            if self._file_format == 'xml':

                utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a SBML model \u27A4\u27A4\u27A4", color="green", style="normal")

                self._biomodel =  self._sbml_reader._read_file(self._file_path)

                if self._biomodel is not None:

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The SBML model: {self._file_name} has been succesfully converted to a BioModel\u27A4\u27A4\u27A4", color="green", style="normal")

            elif self._file_format == 'cellml':

                utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a CellML model \u27A4\u27A4\u27A4", color="green", style="normal")

                self._biomodel = self._cellml_reader._read_file(self._file_path)



    
        

    # ********************************
    # *           Function           *
    # ********************************
    def checkModelReversibility(self, return_irreversibles = False, printing = "off"):

        try:

            if return_irreversibles:

                reversibility, irreversibles = self._model_checker.check_model_reversibility(self._biomodel, return_irreversibles = return_irreversibles)

                if printing.lower() == "on":

                    if reversibility:

                        utility.message_printer("\nAll reactions in the model are REVERSIBLE.")

                    else:

                        utility.message_printer("\nAll reactions in the model are NOT reversible!")

                return reversibility, irreversibles
            
            else:

                reversibility = self._model_checker.check_model_reversibility(self._biomodel)

                if printing.lower() == "on":

                    if reversibility:

                        utility.message_printer("\nAll reactions in the model are REVERSIBLE.")

                    else:

                        utility.message_printer("\nAll reactions in the model are NOT reversible!")

                return reversibility

        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "checkModelConsistency")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")


    # ********************************
    # *           Function           *
    # ********************************
    def getListOfReactions(self):

        return self._biomodel._reactions
    

    # ********************************
    # *           Function           *
    # ********************************
    def getListOfSpecies(self):

        return self._biomodel._species

    
    # ********************************
    # *           Function           *
    # ********************************
    def getStoichiometricMatrix(self, printing = "off"):

        try:
            stoichiometric_matrix = self._matrix_constructor.stoichiometric_matrix_constructor(self._biomodel)
            

            if printing.lower() == "on":
                utility.printer("\nThe Stoichiometric Matrix is:\n", stoichiometric_matrix)
                #time.sleep(2)

            return stoichiometric_matrix

        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return

        except exceptions.EmptyList as e:
            utility.printer("\nAn error has been raised in function: ", "getStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except Exception as e:
            utility.printer("\nAn error has been raised in function: ", "getStoichiometricMatrix")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return


    # ********************************
    # *           Function           *
    # ********************************    
    def getForwardStoichiometricMatrix(self, printing = "off"):

        try:
            forward_stoichiometric_matrix = self._matrix_constructor.forward_stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Forward Stoichiometric Matrix is:\n", forward_stoichiometric_matrix)
                #time.sleep(5)

            return forward_stoichiometric_matrix

        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getForwardStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except exceptions.EmptyList as e:
            utility.printer("\nAn error has been raised in function: ", "getForwardStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            print("Unable to complete the query\!")
            return
        
        except Exception as e:
            utility.printer("\nAn error has been raised in function: ", "getForwardStoichiometricMatrix")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        

    # ********************************
    # *           Function           *
    # ********************************
    def getReverseStoichiometricMatrix(self, printing = "off"):

        try:
            reverse_stoichiometric_matrix = self._matrix_constructor.reverse_stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Reverse Stoichiometric Matrix is:\n", reverse_stoichiometric_matrix)
                #time.sleep(5)

            return reverse_stoichiometric_matrix

        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getReverseStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except exceptions.EmptyList as e:
            utility.printer("\nAn error has been raised in function: ", "getReverseStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except Exception as e:
            utility.printer("\nAn error has been raised in function: ", "getReverseStoichiometricMatrix")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
    
    # ********************************
    # *           Function           *
    # ********************************
    def getStoichiometricColumnNamesIndices(self):

        return self._matrix_constructor.stoichiometric_matrix_column_names(self._biomodel)
    

    # ********************************
    # *           Function           *
    # ********************************
    def getStoichiometricRowNamesIndices(self):

        return self._matrix_constructor.stoichiometric_matrix_row_names(self._biomodel)
    

    # ********************************
    # *           Function           *
    # ********************************
    def getElementInformationInStoichiometricMatrix(self, i, j, printing = "off"):

        try:

            if self._biomodel is None:
                raise exceptions.NoModel("No BioModel has been read!!!")

            return self._matrix_constructor.stoichiometric_matrix_element_information(i, j, self._biomodel, printing = printing)
        
        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getElementInformationInStoichiometricMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        

    # ********************************
    # *           Function           *
    # ********************************
    def getThermoConversionMatrix(self, printing = "off"):

        try:
            conversion_matrix = self._matrix_constructor.kinetic_thermo_conversion_matrix_constructor(self._biomodel, printing = printing)

            return conversion_matrix
            
        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getThermoConversionMatrix")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except Exception as e:
            utility.printer("\nAn error has been raised in function: ", "getThermoConversionMatrix")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return

    
    # ********************************
    # *           Function           *
    # ********************************
    def getKineticRateConstantsVector(self, printing = "off"):

        try:
            kinetic_constants_vector = self._matrix_constructor.kinetic_constants_vector_constructor(self._biomodel, printing)

            return kinetic_constants_vector

        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getKineticRateConstantsVector")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return

        except ValueError as e:
            utility.printer("\nAn error has been raised in function: ", "getKineticRateConstantsVector")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except sp.SympifyError as e:
            utility.printer("\nAn error has been raised in function: ", "getKineticRateConstantsVector")
            utility.error_printer("Sympify Error: ", e)
            utility.message_printer("Equation couldn't be converted to Sympy expression for reaction", color="red", style="normal")
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")

        except exceptions.NoReverseRateConstant as e:
            utility.printer("\nAn error has been raised in function: ", "getKineticRateConstantsVector")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return

        except Exception as e:
            utility.printer("\nAn error has been raised in function: ", "getKineticRateConstantsVector")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")

    
    # ********************************
    # *           Function           *
    # ********************************
    def KineticConstantsThermoCompatibilty(self, printing = "off"):

        try:
            comatibility = self._matrix_constructor.kinetic_rates_thermo_compatibility_check(self._biomodel, printing)

            utility.display_warnings()

            return comatibility
        
        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "KineticConstantsThermoCompatibilty")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return

        except ValueError as e:
            utility.printer("\nAn error has been raised in function: ", "KineticConstantsThermoCompatibilty")
            utility.error_printer("ERROR: ", e)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")
            return
        
        except Exception as e:
            utility.printer("\nAn error has been raised in function: ", "KineticConstantsThermoCompatibilty")
            tb = traceback.extract_tb(sys.exc_info()[2])
            file_path, lineno, _, _ = tb[-1]  # Get the last entry in the traceback
            file_name = os.path.basename(file_path)
            utility.error_printer("Error occurred in file: ", file_name, u_end = "", error_color='yellow')
            utility.error_printer(", line: ", lineno, error_color="yellow")
            utility.error_printer("Unexpected Error: ", e)
            utility.error_printer("Error type: ", sys.exc_info()[0].__name__)
            utility.message_printer("Unable to complete the query\!", color="red", style="normal")