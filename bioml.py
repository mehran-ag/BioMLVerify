import sys
import os
import time

import sbml_reader
import cellml_reader
import matrix_constructor
import model_checker
import exceptions
import utility


from classes.cBioMLReaction import *
from classes.cBioMLModel import *
from classes.cBioMLSpecies import *
from classes.cBioMLParameter import *
from classes.cBioMLSpeciesReference import *
from classes.cBioMLFunctionDefinition import *


class BioML(object):
    '''
    This class creates a model, either CellML or SBML, that will be processed by other classes
    '''

    def __init__(self):
        '''
        Initializes the ModelReader by reading path for a file
        '''
        BioMLSpecies.reset_counter()
        BioMLReaction.reset_counter()
        self._file_path = None
        self._file_name= None
        self._file_format = None
        self._biomlmodel = None
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
    def read_file(self, file_path: str) -> None:
        '''
        Reads the file path and calls smblreader or cellmlreader to convert the input model into a biomlmodel

        Args:
            file_path: A string that defines the path to the file

        Returns:
            None, converts the input model into a BioMLModel, the main internal class, and saves it in an internal variable: self._biomlmodel
        '''

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

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a SBML model \u27A4\u27A4\u27A4\n", color="green")

                    self._biomlmodel =  self._sbml_reader.read_file(self._file_path)

                    file_type = 'SBML'

                elif self._file_format == 'cellml':

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a CellML model \u27A4\u27A4\u27A4\n", color="green")

                    self._biomlmodel = self._cellml_reader.read_file(self._file_path)

                    file_type = 'CellML'

            except Exception as e:
                utility.error_handler(e, function="reading_file")

                return
            
            else:

                if self._biomlmodel is not None:
                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The {file_type} model: {self._file_name} has been succesfully converted to a BioModel \u27A4\u27A4\u27A4\n", color="green")
                else:
                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The imported {file_type} model has not been converted to a BioModel \u27A4\u27A4\u27A4\n", color="red", style="bold")
                    time.sleep(1)







    # ********************************
    # *           Function           *
    # ********************************
    def check_mass_action_kinetics(self, printing = "off"):

        try:

            if self._biomlmodel:

                if self._biomlmodel.is_mass_action is None:

                    is_mass_action = False

                    is_mass_action = self._model_checker.check_mass_action_kinetics(self._biomlmodel)


                    if is_mass_action:

                        if printing.lower() == "on":

                            utility.message_printer(f"\nALL reactions in the model are \"Mass Action\" kinetics\n", color='green')

                            time.sleep(10)

                        return True

                    else:

                        if printing.lower() == "on":

                            utility.message_printer(f"\nModel has (a) reaction(s) not governed by \"Mass Action\" kinetics\n", color='red')        

                            time.sleep(10)

                        return False
                    
                else:

                    if self._biomlmodel.is_mass_action:

                        if printing.lower() == "on":

                            utility.warning_printer(f"\nThis model has been converted from a CellML model and the equations extracted from the model might not be related to equations.\n")

                            utility.message_printer(f"ALL equations in the model are \"Mass Action\" kinetics\n", color='green')

                            time.sleep(5)

                        return True

                    else:

                        utility.warning_printer(f"\nThis model has been converted from a CellML model and the equations extracted from the model might not be related to equations.\n")

                        utility.message_printer(f"\nModel has (a) equation(s) not governed by \"Mass Action\" kinetics\n", color='red')        

                        time.sleep(5)

                        return False

                

            
        except Exception as e:
            utility.error_handler(e, "check_mass_action_kinetics")
            return






    
        

    # ********************************
    # *           Function           *
    # ********************************
    def check_model_reversibility(self, return_irreversibles = False, printing = "off"):

        try:

            if return_irreversibles:

                reversibility, irreversibles = self._model_checker.check_model_reversibility(self._biomlmodel, return_irreversibles = return_irreversibles)

                if printing.lower() == "on":

                    if reversibility is not None:

                        if reversibility:
                            utility.message_printer("\nAll reactions in the model are REVERSIBLE.")
                        else:
                            utility.message_printer("\nAll reactions in the model are NOT reversible!", color='magenta')

                    else:

                        utility.message_printer("\nReversibility conditions are not defined for the reactions!!\n", color='magenta')

                return reversibility, irreversibles
            
            else:

                reversibility = self._model_checker.check_model_reversibility(self._biomlmodel)

                if reversibility is not None:

                    if printing.lower() == "on":
                        if reversibility:
                            utility.message_printer("\nAll reactions in the model are REVERSIBLE.")
                        else:
                            utility.message_printer("\nAll reactions in the model are NOT reversible!")

                else:

                    utility.message_printer("\nReversibility conditions are not defined for the reactions!!\n", color='magenta')

                return reversibility

        except Exception as e:
            utility.error_handler(e, "check_model_reversibility")
            return






    # ********************************
    # *           Function           *
    # ********************************
    def get_list_of_reactions(self):

        return self._biomlmodel._reactions
    






    # ********************************
    # *           Function           *
    # ********************************
    def get_list_of_species(self):

        return self._biomlmodel._species

    




    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_matrix(self, printing = "off"):

        try:
            stoichiometric_matrix = self._matrix_constructor.stoichiometric_matrix_constructor(self._biomlmodel)
            

            if printing.lower() == "on":
                utility.printer("\nThe Stoichiometric Matrix is:\n", stoichiometric_matrix)

            return stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "get_stoichiometric_matrix")
            return






    # ********************************
    # *           Function           *
    # ********************************    
    def get_forward_stoichiometric_matrix(self, printing = "off"):

        try:
            forward_stoichiometric_matrix = self._matrix_constructor.forward_stoichiometric_matrix_constructor(self._biomlmodel)

            if printing.lower() == "on":
                utility.printer("\nThe Forward Stoichiometric Matrix is:\n", forward_stoichiometric_matrix)

            return forward_stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "get_forward_stoichiometric_matrix")
            return
        





    # ********************************
    # *           Function           *
    # ********************************
    def get_reverse_stoichiometric_matrix(self, printing = "off"):

        try:
            reverse_stoichiometric_matrix = self._matrix_constructor.reverse_stoichiometric_matrix_constructor(self._biomlmodel)

            if printing.lower() == "on":
                utility.printer("\nThe Reverse Stoichiometric Matrix is:\n", reverse_stoichiometric_matrix)

            return reverse_stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "get_reverse_stoichiometric_matrix")
            return
        
    




    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_column_names_indices(self, printing="off"):

        columns = self._matrix_constructor.stoichiometric_matrix_column_names(self._biomlmodel)

        if printing.lower() == "on":

            utility.printer("\nThe Columns are:", columns )

        return columns
    






    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_row_names_indices(self, printing="off"):

        rows = self._matrix_constructor.stoichiometric_matrix_row_names(self._biomlmodel)

        if printing.lower() == "on":

            utility.printer("\nThe Rows are:", rows )

        return rows
    






    # ********************************
    # *           Function           *
    # ********************************
    def get_element_information_in_stoichiometric_matrix(self, i, j, printing = "off"):

        try:

            if self._biomlmodel is None:
                raise exceptions.NoModel("No BioModel has been read!!!")

            return self._matrix_constructor.stoichiometric_matrix_element_information(i, j, self._biomlmodel, printing = printing)
        
        except Exception as e:
            utility.error_handler(e, "get_element_information_in_stoichiometric_matrix")
            return
        





    # ********************************
    # *           Function           *
    # ********************************
    def get_thermo_conversion_matrix(self, printing = "off"):

        try:
            conversion_matrix = self._matrix_constructor.kinetic_thermo_conversion_matrix_constructor(self._biomlmodel, printing = printing)

            return conversion_matrix
            
        except Exception as e:
            utility.error_handler(e, "get_thermo_conversion_matrix")
            return

    




    # ********************************
    # *           Function           *
    # ********************************
    def get_kinetic_rate_constants_vector(self, printing = "off"):

        try:
            kinetic_constants_vector = self._matrix_constructor.kinetic_constants_vector_constructor(self._biomlmodel, printing)

            return kinetic_constants_vector
        
        except Exception as e:
            utility.error_handler(e, "get_kinetic_rate_constants_vector")
            return

    



    # ********************************
    # *           Function           *
    # ********************************
    def check_kinetic_constants_thermo_compatibility(self, printing = "off"):

        try:
            compatibility = self._matrix_constructor.kinetic_rates_thermo_compatibility_check(self._biomlmodel, printing)

            utility.display_warnings()

            return compatibility
        
        except Exception as e:
            utility.error_handler(e, "check_kinetic_constants_thermo_compatibility")
            utility.printer("\nCompatibility Check: ","\nThe kinetic reaction rate constants are NOT compatible with thermodynamic constraints\n", text_color="red", text_style="bold")
            return