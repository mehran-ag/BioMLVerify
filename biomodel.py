import sys
import os
import time

import sbml_reader
import cellml_reader
import matrix_constructor
import model_checker
import exceptions
import utility


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
        Initializes the ModelReader by reading path for a file
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







    # ********************************
    # *           Function           *
    # ********************************
    def read_file(self, file_path):

        try:
            if not isinstance(file_path, str):
                raise TypeError("File path must be a string.")
            
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
            
            self._file_path = file_path
            self._file_name = os.path.basename(file_path)
            self._file_format = os.path.splitext(file_path)[1][1:]
        
        except Exception as e:
            utility.error_handler(e)
            return

        else:

            try:
                if self._file_format == 'xml':

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a SBML model \u27A4\u27A4\u27A4", color="green", style="normal")

                    self._biomodel =  self._sbml_reader.read_file(self._file_path)

                    file_type = 'SBML'

                elif self._file_format == 'cellml':

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a CellML model \u27A4\u27A4\u27A4", color="green", style="normal")

                    self._biomodel = self._cellml_reader.read_file(self._file_path)

                    file_type = 'CellML'

            except Exception as e:
                utility.error_handler(e, function="reading_file")
                sys.exit("\n\nExecution terminated as there is no file read to continue process!!\n\n")
            
            else:

                if self._biomodel is not None:
                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The {file_type} model: {self._file_name} has been succesfully converted to a BioModel\u27A4\u27A4\u27A4", color="green", style="normal")
                else:
                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The imported {file_type} model has not been converted to a BioModel\u27A4\u27A4\u27A4", color="red", style="bold")
                    time.sleep(3)







    # ********************************
    # *           Function           *
    # ********************************
    def checkMassActionKinetics(self, printing = "off"):

        try:

            is_mass_action = self._model_checker.check_mass_action_kinetics(self._biomodel)

            if is_mass_action:

                if printing.lower() == "on":

                    utility.message_printer(f"ALL reactions in the model are \"Mass Action\" kinetics\n\n\n", color='green', style='normal')

                    time.sleep(10)

                return True

            else:

                if printing.lower() == "on":

                    utility.message_printer(f"Model has (a) reaction(s) not governed by \"Mass Action\" kinetics\n\n\n", color='red', style='normal')        

                    time.sleep(10)

                return False
            
        except Exception as e:
            utility.error_handler(e, "checkMassActionKinetics")
            return






    
        

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

        except Exception as e:
            utility.error_handler(e, "checkModelReversibility")
            return






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

            return stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "getStoichiometricMatrix")
            return






    # ********************************
    # *           Function           *
    # ********************************    
    def getForwardStoichiometricMatrix(self, printing = "off"):

        try:
            forward_stoichiometric_matrix = self._matrix_constructor.forward_stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Forward Stoichiometric Matrix is:\n", forward_stoichiometric_matrix)

            return forward_stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "getForwardStoichiometricMatrix")
            return
        





    # ********************************
    # *           Function           *
    # ********************************
    def getReverseStoichiometricMatrix(self, printing = "off"):

        try:
            reverse_stoichiometric_matrix = self._matrix_constructor.reverse_stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Reverse Stoichiometric Matrix is:\n", reverse_stoichiometric_matrix)

            return reverse_stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "getReverseStoichiometricMatrix")
            return
        
    




    # ********************************
    # *           Function           *
    # ********************************
    def getStoichiometricColumnNamesIndices(self, printing="off"):

        columns = self._matrix_constructor.stoichiometric_matrix_column_names(self._biomodel)

        if printing.lower() == "on":

            utility.printer("\nThe Columns are:", columns )

        return columns
    






    # ********************************
    # *           Function           *
    # ********************************
    def getStoichiometricRowNamesIndices(self, printing="off"):

        rows = self._matrix_constructor.stoichiometric_matrix_row_names(self._biomodel)

        if printing.lower() == "on":

            utility.printer("\nThe Rows are:", rows )

        return rows
    






    # ********************************
    # *           Function           *
    # ********************************
    def getElementInformationInStoichiometricMatrix(self, i, j, printing = "off"):

        try:

            if self._biomodel is None:
                raise exceptions.NoModel("No BioModel has been read!!!")

            return self._matrix_constructor.stoichiometric_matrix_element_information(i, j, self._biomodel, printing = printing)
        
        except Exception as e:
            utility.error_handler(e, "getElementInformationInStoichiometricMatrix")
            return
        





    # ********************************
    # *           Function           *
    # ********************************
    def getThermoConversionMatrix(self, printing = "off"):

        try:
            conversion_matrix = self._matrix_constructor.kinetic_thermo_conversion_matrix_constructor(self._biomodel, printing = printing)

            return conversion_matrix
            
        except Exception as e:
            utility.error_handler(e, "getThermoConversionMatrix")
            return

    




    # ********************************
    # *           Function           *
    # ********************************
    def getKineticRateConstantsVector(self, printing = "off"):

        try:
            kinetic_constants_vector = self._matrix_constructor.kinetic_constants_vector_constructor(self._biomodel, printing)

            return kinetic_constants_vector
        
        except Exception as e:
            utility.error_handler(e, "getKineticRateConstantsVector")
            return

    



    # ********************************
    # *           Function           *
    # ********************************
    def KineticConstantsThermoCompatibilty(self, printing = "off"):

        try:
            comatibility = self._matrix_constructor.kinetic_rates_thermo_compatibility_check(self._biomodel, printing)

            utility.display_warnings()

            return comatibility
        
        except Exception as e:
            utility.error_handler(e, "KineticConstantsThermoCompatibilty")
            utility.printer("\nCompatibility Check: ","\nThe kinetic reaction rate constants are NOT compatible with thermodynamic constraints\n", text_color="red", text_style="bold")
            return