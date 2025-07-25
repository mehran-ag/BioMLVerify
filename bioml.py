import os
import time

import sbml_reader
import cellml_reader
import matrix_constructor
import model_checker
import exceptions
import utility

import numpy as np


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
                file_path (str): A string that defines the path to the file

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

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a SBML model \u27A4\u27A4\u27A4\n", color="cyan")

                    self._biomlmodel =  self._sbml_reader.read_file(self._file_path)

                    file_type = 'SBML'

                elif self._file_format == 'cellml':

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a CellML model \u27A4\u27A4\u27A4\n", color="cyan")

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
    def check_mass_action_kinetics(self, printing: bool = False) -> bool:
        '''
            Checks the equations of the imported model to confirm if all equations in the model are governed by Mass Action Kinetics

            Args:
                None

            Returns:
                bool: Boolean value showing if all reactions are governed by Mass Action Kinetics or not
                None: If an exception is raised during the check.
        '''

        try:

            if self._biomlmodel:

                if self._biomlmodel.is_mass_action is None:

                    is_mass_action = False

                    is_mass_action = self._model_checker.check_mass_action_kinetics(self._biomlmodel)


                    if is_mass_action:

                        if printing:

                            utility.message_printer(f"\nALL reactions in the model are \"Mass Action\" kinetics\n", color='green')

                            time.sleep(5)

                        return True

                    else:

                        if printing:

                            utility.message_printer(f"\nModel has (a) reaction(s) not governed by \"Mass Action\" kinetics\n", color='red')        

                            time.sleep(5)

                        return False
                    
                else:

                    if self._biomlmodel.is_mass_action:

                        if printing:

                            utility.message_printer(f"ALL equations in the model are \"Mass Action\" kinetics\n", color='green')

                            time.sleep(5)

                        return True

                    else:

                        if printing:

                            utility.message_printer(f"\nModel has (a) equation(s) not governed by \"Mass Action\" kinetics\n", color='red')        

                            time.sleep(5)

                        return False
                    
            else:
                
                raise ValueError("No BioMLModel provided")

                

            
        except Exception as e:
            utility.error_handler(e, "check_mass_action_kinetics")
            return None






    
        

    # ********************************
    # *           Function           *
    # ********************************
    def check_model_reversibility(self, return_irreversibles: bool = False, printing: bool = False) -> bool:
        '''
            Checks all reactions of the imported model to confirm if the reactions are Reversibe

            Args:
                return_irreversibles (bool): if this value is True, a list containing the irreversible reactions will be returned
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                bool: Boolean value showing if all reactions are reversible or not
                None: If an exception is raised during the check.
        '''

        try:

            if return_irreversibles:

                reversibility, irreversibles = self._model_checker.check_model_reversibility(self._biomlmodel, return_irreversibles = return_irreversibles)

                if printing:

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

                    if printing:
                        if reversibility:
                            utility.message_printer("\nAll reactions in the model are REVERSIBLE.")
                        else:
                            utility.message_printer("\nAll reactions in the model are NOT reversible!")

                else:

                    utility.message_printer("\nReversibility conditions are not defined for the reactions!!\n", color='magenta')

                return reversibility

        except Exception as e:
            utility.error_handler(e, "check_model_reversibility")
            return None






    # ********************************
    # *           Function           *
    # ********************************
    def get_list_of_reactions(self) -> list[object]:
        '''
            Returns the list of BioML Reactions (BioML reaction instances) in the model

            Args:
                None

            Returns:
                list: A list containing the BioML reaction instances
        '''

        return self._biomlmodel._reactions
    






    # ********************************
    # *           Function           *
    # ********************************
    def get_list_of_species(self) -> list[object]:
        '''
            Returns the list of BioMLSpecies (BioML species instances) in the model

            Args:
                None

            Returns:
            list: A list containing the BioML species instances
        '''

        return self._biomlmodel._species

    




    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_matrix(self, printing: bool = False) -> np.ndarray:
        '''
            Returns a 2D numpy array representing the stoichiometric matrix of the reactions for the imported model: This matrix represents the stoichiometric coefficients of species in reactions having species as rows and reactions as columns

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A numpy array
                None: If an exception is raised during the execution.
        '''

        try:
            stoichiometric_matrix = self._matrix_constructor.construct_stoichiometric_matrix(self._biomlmodel)
            

            if printing:
                utility.printer("\nThe Stoichiometric Matrix is:\n", stoichiometric_matrix)

            return stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "get_stoichiometric_matrix")
            return None






    # ********************************
    # *           Function           *
    # ********************************    
    def get_forward_stoichiometric_matrix(self, printing: bool = False) -> np.ndarray:
        '''
            Returns a 2D numpy array representing the forward stoichiometric matrix of the reactions for the imported model: This array shows stoichiometric coefficients of reactants only (all values are positive in this array)

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 2D numpy array
                None: If an exception is raised during the execution.
        '''

        try:
            forward_stoichiometric_matrix = self._matrix_constructor.construct_forward_stoichiometric_matrix(self._biomlmodel)

            if printing:
                utility.printer("\nThe Forward Stoichiometric Matrix is:\n", forward_stoichiometric_matrix)

            return forward_stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "get_forward_stoichiometric_matrix")
            return None
        





    # ********************************
    # *           Function           *
    # ********************************
    def get_reverse_stoichiometric_matrix(self, printing: bool = False) -> np.ndarray:
        '''
            Returns a 2D numpy array representing the reverse stoichiometric matrix of the reactions for the imported model: This array shows stoichiometric coefficients of products only (all values are positive in this array)

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 2D numpy array
                None: If an exception is raised during the execution.
        '''

        try:
            reverse_stoichiometric_matrix = self._matrix_constructor.construct_reverse_stoichiometric_matrix(self._biomlmodel)

            if printing:
                utility.printer("\nThe Reverse Stoichiometric Matrix is:\n", reverse_stoichiometric_matrix)

            return reverse_stoichiometric_matrix

        except Exception as e:
            utility.error_handler(e, "get_reverse_stoichiometric_matrix")
            return None
        
    




    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_column_names_indices(self, printing: bool = False) -> dict:
        '''
            Returns a dictionary containing the names of columns with their corresponding indices in the stoichiometric matrix

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                dict: A dictionary mapping column index to column name
        '''

        columns = self._matrix_constructor.get_stoichiometric_matrix_column_names(self._biomlmodel)

        if printing:

            utility.printer("\nThe Columns are:", columns )

        return columns
    






    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_row_names_indices(self, printing: bool = False) -> dict:
        '''
            Returns a dictionary containing the names of rows with their corresponding indices in the stoichiometric matrix

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                dict: A dictionary mapping row index to row name
        '''

        rows = self._matrix_constructor.get_stoichiometric_matrix_row_names(self._biomlmodel)

        if printing:

            utility.printer("\nThe Rows are:", rows )

        return rows
    






    # ********************************
    # *           Function           *
    # ********************************
    def get_element_information_in_stoichiometric_matrix(self, i: int, j: int, printing: bool = False) -> str:
        '''
            Returns the element in the stoichiometric matrix

            Args:
                i (int): the row index
                j (int): the column index 
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                str: A value (string)
                None: If an exception is raised during the execution.
        '''

        try:

            if self._biomlmodel is None:
                raise exceptions.NoModel("No BioModel has been read!!!")

            return self._matrix_constructor.get_stoichiometric_matrix_element_information(i, j, self._biomlmodel, printing = printing)
        
        except Exception as e:
            utility.error_handler(e, "get_element_information_in_stoichiometric_matrix")
            return None
        





    # ********************************
    # *           Function           *
    # ********************************
    def get_thermo_conversion_matrix(self, printing: bool = False) -> np.ndarray:
        '''
            Returns a 2D numpy array that converts kinetic reaction rate constants to corresponding thermodynamic reaction rate constants

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 2D numpy array
                None: If an exception is raised during the execution.
        '''

        try:
            conversion_matrix = self._matrix_constructor.construct_kinetic_thermo_conversion_matrix(self._biomlmodel, printing = printing)

            return conversion_matrix
            
        except Exception as e:
            utility.error_handler(e, "get_thermo_conversion_matrix")

            return None

    




    # ********************************
    # *           Function           *
    # ********************************
    def get_kinetic_rate_constants_vector(self, printing: bool = False) -> np.ndarray:
        '''
            Returns a 1D numpy array that contains the ratio of forward to reverse reaction rate constants

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 1D numpy array
                None: If an exception is raised during the execution.
        '''

        try:
            kinetic_constants_vector = self._matrix_constructor.construct_kinetic_constants_vector(self._biomlmodel, printing)

            return kinetic_constants_vector
        
        except Exception as e:
            utility.error_handler(e, "get_kinetic_rate_constants_vector")

            return None

    



    # ********************************
    # *           Function           *
    # ********************************
    def check_kinetic_constants_thermo_compatibility(self, printing: bool = False) -> bool:
        '''
            Checks the validity of Kinetic reaction rate constants in thermodynamic framework.
            The function uses Wegscheider conditions to check thermodynamic compatibility of constants

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                bool: True if the reaction rate constants are compatible and meaningful, False otherwise
                None: If an exception is raised during the execution.
        '''

        try:
            compatibility = self._matrix_constructor.check_kinetic_rates_thermo_compatibility(self._biomlmodel, printing)

            utility.display_warnings()

            return compatibility
        
        except Exception as e:
            utility.error_handler(e, "check_kinetic_constants_thermo_compatibility")
            utility.printer("\nThermodynamic Compatibility Check: ","\nThe kinetic reaction rate constants are NOT compatible with thermodynamic constraints\n", text_color="red", text_style="bold")
            return None
        


    # ********************************
    # *           Function           *
    # ********************************
    def get_elemental_matrix(self, printing: bool = False) -> np.ndarray:
        '''
            Returns a 2D numpy array that represents the composition of each compound in an array: rows as chemical elements and columns as chemical compounds

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 2D numpy array
                None: If an exception is raised during the execution.
        '''

        try:
            elemental_matrix = self._matrix_constructor.construct_elemental_matrix(self._biomlmodel)
            

            if printing:
                utility.printer("\nThe Stoichiometric Matrix is:\n", elemental_matrix)

            return elemental_matrix

        except Exception as e:
            utility.error_handler(e, "get_elemental_matrix")
            return None
        


    # ********************************
    # *           Function           *
    # ********************************
    def check_mass_balance(self, printing: bool = False) -> bool:
        '''
            Returns True if mass is conserved in all reactions of the model and False if not conserved. If printing is on, it can display the reaction violating mass conservation if mass balance fails

            Args:
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                bool: True if mass is conserved, False otherwise
                None: If an exception is raised during the check.
        '''


        try:

            elemental_array = self.get_elemental_matrix()

            if elemental_array is None:
                raise ValueError(f"Elemental matrix is not provided!")

            stoichiometric_array = self.get_stoichiometric_matrix()

            if stoichiometric_array is None:
                raise ValueError(f"Stoichiometric matrix is not provided!")

            if elemental_array.shape[1] == stoichiometric_array.shape[0]:
                conservation_array = elemental_array @ stoichiometric_array
            
            else:

                raise ValueError(f"Matrices cannot be multiplied as the number of columns of the elemental matrix, {elemental_array.shape[1]}, is not equal to the number of rows of the stoichiometric matrix, {stoichiometric_array.shape[0]}!")
            
            if np.all( ( conservation_array == 0 ) ):

                if printing:
                    utility.message_printer("\nMass is conserved in the reactions\n", color='green')

                return True
        
            else:
                if printing:
                    utility.message_printer("\nConservation of Mass is violated", color='red')
                
                    length = conservation_array.shape[1]
                    for i in  range(0,length):
                        if np.all(conservation_array[:, i] == 0):
                            pass
                        else:
                            for biomlreaction in self._biomlmodel.reactions:
                                if biomlreaction.index == i:
                                    utility.message_printer(f"\nMass is not conserved in reaction {biomlreaction.ID}", color='magenta')

                return False


        
        except Exception as e:
            utility.error_handler(e, "check_mass_balance")
            return None
        



    def verify_model(self, mass_balance: bool = False, printing: bool = False):
        '''
        Returns True if model complies with thermodynamic principles and False if not.

        Args:
            mass_balance (bool): if this value is True, mass conservation in reactions can also be checked
            printing: if this value is True, a message will be displayed to show the result

        Returns:
           bool: True if model is consistent with thermodynamic rules, False otherwise
           None: If an exception is raised during the check.
        '''

        passed = True

        try:

            is_mass_balanced = None

            if mass_balance:

                is_mass_balanced = self.check_mass_balance(printing)

                if is_mass_balanced is None:
                    raise ValueError(f"Mass balance cannot be checked!")

                if is_mass_balanced:

                    if printing:

                        utility.printer("\nMass Conservation Check: ",f"Mass is conserved in the reactions\n", text_color="green", text_style="bold")

                else:

                    passed = False
                        
                    if printing:

                        utility.printer("\nMass Conservation Check: ",f"Mass is NOT conserved in the reactions", text_color="red", text_style="bold")

        except Exception as e:
            utility.error_handler(e, "verify_model")


        try:

            if passed:

                if not self.check_mass_action_kinetics():

                    passed = False

                    if printing:

                        utility.printer("\nThermodynamic Compatibility Check: ",f"Model {self._file_name} has (a) reaction(s) not governed by \"Mass Action\" kinetics and\n{' ' * 37}it is NOT eligible for verification check\n", text_color="red")

                    return passed
                
                if not self.check_model_reversibility():

                    passed = False

                    if printing:

                        utility.printer("\nThermodynamic Compatibility Check: ",f"Model {self._file_name} has (an) irrversible reaction(s) and\n{' ' * 37}it is NOT eligible for verification check\n", text_color="red")

                    return passed
                
                if self.check_kinetic_constants_thermo_compatibility(printing):
                    passed = True

                    return passed

        except Exception as e:
            utility.error_handler(e, "verify_model")
            return None