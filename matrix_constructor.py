import numpy as np
import exceptions
from classes.cBioMLModel import BioMLModel
from classes.cBioMLReaction import *
from classes.cBioMLSpecies import *
import utility
import warnings
from constants import *

from scipy.linalg import null_space
    



class MatrixConstructor:


    # ********************************
    # *           Function           *
    # ********************************
    def construct_stoichiometric_matrix(self, biomlmodel: BioMLModel) -> np.ndarray:
        """
            Constructs the stoichiometric matrix for the given BioModel.

            Parameters:
                biomlmodel (BioMLModel): A model of BioML class containing species and reactions.

            Returns:
                np.ndarray: A 2D array representing the stoichiometric matrix, 
                            where rows correspond to species and columns to reactions.
        """

        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomlmodel.get_list_of_species()
        reactions_list = biomlmodel.get_list_of_reactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = BioMLSpecies.get_current_index()

        columns = BioMLReaction.get_current_index()

        self.stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_reactants = individual_reaction.get_list_of_reactants()

            for individual_reactant in reaction_reactants:

                row = individual_reactant.index

                stoichiometry = individual_reactant.get_stoichiometry()

                self.stoichiometric_matrix[row, column] = -1 * int(stoichiometry)

            reaction_products = individual_reaction.get_list_of_products()

            for individual_product in reaction_products:

                row = individual_product.index

                stoichiometry = individual_product.get_stoichiometry()

                self.stoichiometric_matrix[row, column] = int(stoichiometry)

        return self.stoichiometric_matrix
    


    # ********************************
    # *           Function           *
    # ********************************
    def construct_forward_stoichiometric_matrix(self, biomlmodel: BioMLModel) -> np.ndarray:
        """
            Constructs the forward stoichiometric matrix for the given BioModel.

            Parameters:
                biomlmodel (an instance of BioModel class): The biological model containing species and reactions.

            Returns:
                np.ndarray: A 2D array representing the stoichiometric matrix, 
                            where rows correspond to species and columns to reactions.
        """

        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomlmodel.get_list_of_species()
        parameters_list = biomlmodel.get_list_of_parameters()
        reactions_list = biomlmodel.get_list_of_reactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = BioMLSpecies.get_current_index()

        columns = BioMLReaction.get_current_index()

        self.forward_stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_reactants = individual_reaction.get_list_of_reactants()

            for individual_reactant in reaction_reactants:

                row = individual_reactant.index

                stoichiometry = individual_reactant.get_stoichiometry()

                self.forward_stoichiometric_matrix[row, column] = int(stoichiometry)

        return self.forward_stoichiometric_matrix
    


    # ********************************
    # *           Function           *
    # ********************************
    def construct_reverse_stoichiometric_matrix(self, biomlmodel: BioMLModel) -> np.ndarray:
        """
            Constructs the reverse stoichiometric matrix for the given BioModel.

            Parameters:
                biomlmodel (an instance of BioModel class): The biological model containing species and reactions.

            Returns:
                np.ndarray: A 2D array representing the stoichiometric matrix, 
                            where rows correspond to species and columns to reactions.
        """

        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomlmodel.get_list_of_species()
        parameters_list = biomlmodel.get_list_of_parameters()
        reactions_list = biomlmodel.get_list_of_reactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = BioMLSpecies.get_current_index()

        columns = BioMLReaction.get_current_index()

        self.reverse_stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_products = individual_reaction.get_list_of_products()

            for individual_product in reaction_products:

                row = individual_product.index

                stoichiometry = individual_product.get_stoichiometry()

                self.reverse_stoichiometric_matrix[row, column] = int(stoichiometry)

        return self.reverse_stoichiometric_matrix


    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_matrix_column_names(self, biomlmodel: BioMLModel) -> dict:
        """
            Returns a dictionary containing the names of columns with their corresponding indices in the stoichiometric matrix

            Args:
                biomlmodel (BioMlModel): A model of BioML class containing species and reactions.

            Returns:
                dict: A dictionary mapping column index to column name
        """
        
        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        
        reactions_list = biomlmodel.get_list_of_reactions()

        column_indices_names = {}

        for reaction in reactions_list:

            column_indices_names[reaction.index] = reaction.ID

        return column_indices_names
    



    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_matrix_row_names(self, biomlmodel: BioMLModel) -> dict:
        """
            Returns a dictionary containing the names of rows with their corresponding indices in the stoichiometric matrix

            Args:
                biomlmodel (BioMlModel): A model of BioML class containing species and reactions.

            Returns:
                dict: A dictionary mapping row index to row name
        """
    
        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomlmodel.get_list_of_species()

        row_indices_names = {}

        for species in species_list:

            if species.compound:

                row_indices_names[species.index] = species.compound

            else:

                row_indices_names[species.index] = species.ID

        return row_indices_names        
    



    # ********************************
    # *           Function           *
    # ********************************
    def get_stoichiometric_matrix_element_information(self, i: int, j: int, biomlmodel: BioMLModel, printing: bool = False) -> str:
        """
            Returns the element in the stoichiometric matrix

            Args:
                i (int): the row index
                j (int): the column index
                biomlmodel (BioMlModel): A model of BioML class containing species and reactions.
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                str: A value (string)
        """

        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        highest_i = self.stoichiometric_matrix.shape[0]
        highest_j = self.stoichiometric_matrix.shape[1]

        if i > highest_i or j > highest_j:

            utility.error_printer("Invlaid indices provided!")
            
            print(f"\nThe highest value for i (rows), j (columns) are {highest_i} and {highest_j}, respectively")
            return

        else:

            row_indices_names = self.get_stoichiometric_matrix_row_names(biomlmodel)
            species = row_indices_names[i]

            column_indices_names = self.get_stoichiometric_matrix_column_names(biomlmodel)
            reaction = column_indices_names[j]

            if printing:
                utility.printer(f"\nThe stoichiometric coefficient for {species} in reaction {reaction} is: ", f"{self.stoichiometric_matrix[i][j]}")

            return f"The stoichiometric coefficient for {species} in reaction {reaction} is: {self.stoichiometric_matrix[i][j]}"
    


    # ********************************
    # *           Function           *
    # ********************************
    def construct_kinetic_constants_vector(self, biomlmodel: BioMLModel, printing: bool = False) -> np.ndarray:
        """
            Returns a 1D numpy array that contains the ratio of forward to reverse reaction rate constants

            Args:
            biomlmodel (BioMlModel): A model of BioML class containing species and reactions.
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 1D numpy array
        """

        if biomlmodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        
        biomlmodel_reactions = biomlmodel.get_list_of_reactions()

        reactions_number = BioMLReaction.get_current_index()

        kinetic_vector_length = reactions_number

        vector_of_kinetic_constants = np.zeros(kinetic_vector_length)

        for biomlmodel_reaction in biomlmodel_reactions:

            if biomlmodel_reaction.index != None:

                index = biomlmodel_reaction.index

                name = biomlmodel_reaction.get_id()

                k_plus_value = biomlmodel_reaction.kinetic_forward_rate_constant_value if biomlmodel_reaction.kinetic_forward_rate_constant_value is not None else 0.0

                k_minus_value = biomlmodel_reaction.kinetic_reverse_rate_constant_value if biomlmodel_reaction.kinetic_reverse_rate_constant_value is not None else 0.0

                if k_minus_value == 0.:
                    if biomlmodel_reaction.kinetic_reverse_rate_constant:
                        raise exceptions.NoReverseRateConstant(f"Kinetic Constants Vector cannot be constructed since there is no initial value (or it is zero) for the reverse reaction rate constant for reaction {name}: {biomlmodel_reaction.get_kinetic_law()}")
                    else:
                        raise exceptions.NoReverseRateConstant(f"Kinetic Constants Vector cannot be constructed since there is no reverse reaction rate constant for reaction {name}: {biomlmodel_reaction.get_kinetic_law()}")

                vector_of_kinetic_constants[index] = k_plus_value / k_minus_value

        if printing:
            utility.printer("\nKinetic Constants Vector is:\n",vector_of_kinetic_constants)

        return vector_of_kinetic_constants
    

    # ********************************
    # *           Function           *
    # ********************************
    def construct_kinetic_thermo_conversion_matrix(self, biomlmodel: BioMLModel, printing: bool = False) -> np.ndarray:
        """
            Returns a 2D numpy array that converts kinetic reaction rate constants to corresponding thermodynamic reaction rate constants

            Args:
                biomlmodel (BioMlModel): A model of BioML class containing species and reactions.
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                np.ndarray: A 2D numpy array
        """

        if biomlmodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        reactions_number = BioMLReaction.get_current_index()

        identity_array = np.eye(reactions_number)

        forward_stoichiometric_matrix = self.construct_forward_stoichiometric_matrix(biomlmodel)

        reverse_stoichiometric_matrix = self.construct_reverse_stoichiometric_matrix(biomlmodel)

        transposed_forward_stoichiometric_matrix = np.transpose(forward_stoichiometric_matrix)

        transposed_reverse_stoichiometric_matrix = np.transpose(reverse_stoichiometric_matrix)

        conversion_matrix = np.block( [ [ identity_array, transposed_forward_stoichiometric_matrix ], [ identity_array, transposed_reverse_stoichiometric_matrix ] ] )

        if printing:
            utility.printer("\nConversion Matrix is\n:", conversion_matrix)

        return conversion_matrix
    

    # ********************************
    # *           Function           *
    # ********************************
    def check_kinetic_rates_thermo_compatibility(self, biomlmodel: BioMLModel, printing: bool = False) -> bool:
        """
            Checks the validity of Kinetic reaction rate constants in thermodynamic framework.
            The function uses Wegscheider conditions to check thermodynamic compatibility of constants

            Args:
                biomlmodel (BioMlModel): A model of BioML class containing species and reactions.
                printing (bool): if this value is True, a message will be displayed to show the result

            Returns:
                bool: True if the reaction rate constants are compatible and meaningful, False otherwise
        """

        if biomlmodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        

        kinetic_rates_vector = self.construct_kinetic_constants_vector(biomlmodel)

        with warnings.catch_warnings():

            warnings.simplefilter('error', RuntimeWarning)

            logn_kinetic_rates_vector = np.log(kinetic_rates_vector)

        logn_kinetic_rates_vector = logn_kinetic_rates_vector.reshape(-1,1)

        stoichiometric_matrix = self.construct_stoichiometric_matrix(biomlmodel)

        minus_stoichiometric_matrix = -1 * stoichiometric_matrix

        minus_stoichiometric_null_space = null_space(minus_stoichiometric_matrix)

        transposed_minus_stoichiometric_null_space = minus_stoichiometric_null_space.T

        result = transposed_minus_stoichiometric_null_space @ logn_kinetic_rates_vector

        if np.all(np.abs(result) <= 1e-2):

            if printing:
                utility.printer("\nThermodynamic Compatibility Check: ","The kinetic reaction rate constants are compatible with thermodynamic constraints\n", text_color="green", text_style="bold")

            return True
        
        else:

            if printing:
                utility.printer("\nThermodynamic Compatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints\n", text_color="red", text_style="bold")
            
            return False
        



    # ********************************
    # *           Function           *
    # ********************************
    def construct_elemental_matrix(self, biomlmodel: BioMLModel) -> np.ndarray:
        """
        Constructs the elemental matrix for the given BioModel.

        Parameters:
            biomlmodel (BioMlModel): A model of BioML class containing species and reactions.

        Returns:
            np.ndarray: A 2D array representing the elemental matrix, 
                        where rows correspond to elements and columns to species.
        """
        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        biomlspecies_list = biomlmodel.get_list_of_species()

        if len(biomlspecies_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        

        columns = BioMLSpecies.get_current_index()

        element_indices_dict = biomlmodel.mk_element_indices_dict()

        rows = len( element_indices_dict )

        self.elemental_matrix = np.zeros((rows, columns), dtype = int)

        for individual_biomlspecies in biomlspecies_list:

            column = individual_biomlspecies.index

            for compound, number in individual_biomlspecies.composition.items():

                row = element_indices_dict.get(compound)

                if row is None:
                    raise ValueError(f"There is not an index for {compound} in species {individual_biomlspecies.name if individual_biomlspecies.name is not None else individual_biomlspecies.ID}")
                
                self.elemental_matrix[row, column] = number
                

        return self.elemental_matrix
    





    # ********************************
    # *           Function           *
    # ********************************
    def construct_charge_matrix(self, biomlmodel: BioMLModel) -> np.ndarray:
        """
            Constructs the charge matrix for the given BioModel.

            Parameters:
                biomlmodel (BioMLModel): A model of BioML class containing species and reactions.

            Returns:
                np.ndarray: A 2D array representing the stoichiometric matrix, 
                            where rows correspond to species and columns to reactions.
        """

        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomlmodel.get_list_of_species()
        reactions_list = biomlmodel.get_list_of_reactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = BioMLSpecies.get_current_index()

        columns = BioMLReaction.get_current_index()

        charge_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_reactants = individual_reaction.get_list_of_reactants()

            for individual_reactant in reaction_reactants:

                row = individual_reactant.index

                charge = individual_reactant.get_charge()

                charge_matrix[row, column] = charge

            reaction_products = individual_reaction.get_list_of_products()

            for individual_product in reaction_products:

                row = individual_product.index

                charge = individual_product.get_charge()

                charge_matrix[row, column] = int(charge)

        transposed_charge_matrix = np.transpose(charge_matrix)

        return transposed_charge_matrix


