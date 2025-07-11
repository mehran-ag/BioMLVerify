import numpy as np
import exceptions
import os
import sympy as sp
from classes.cReaction import *
from classes.cSpecies import *
import utility
import random
import warnings
import re
import inspect
from constants import *
from collections import defaultdict

from scipy.linalg import null_space
    



class MatrixConstructor:


    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_constructor(self, biomodel) -> np.ndarray:
        """
        Constructs the stoichiometric matrix for the given BioModel.

        Parameters:
            biomodel (an instance of BioModel class): The biological model containing species and reactions.

        Returns:
            np.ndarray: A 2D array representing the stoichiometric matrix, 
                        where rows correspond to species and columns to reactions.
        """

        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomodel.getListOfSpecies()
        parameters_list = biomodel.getListOfParameters()
        reactions_list = biomodel.getListOfReactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = Species.getCurrentIndex()

        columns = Reaction.getCurrentIndex()

        self.stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_reactants = individual_reaction.getListOfReactants()

            for individual_reactant in reaction_reactants:

                row = individual_reactant.index

                stoichiometry = individual_reactant.getStoichiometry()

                self.stoichiometric_matrix[row, column] = -1 * int(stoichiometry)

            reaction_products = individual_reaction.getListOfProducts()

            for individual_product in reaction_products:

                row = individual_product.index

                stoichiometry = individual_product.getStoichiometry()

                self.stoichiometric_matrix[row, column] = int(stoichiometry)

        return self.stoichiometric_matrix
    


    # ********************************
    # *           Function           *
    # ********************************
    def forward_stoichiometric_matrix_constructor(self, biomodel) -> np.ndarray:
        """
        Constructs the forward stoichiometric matrix for the given BioModel.

        Parameters:
            biomodel (an instance of BioModel class): The biological model containing species and reactions.

        Returns:
            np.ndarray: A 2D array representing the stoichiometric matrix, 
                        where rows correspond to species and columns to reactions.
        """

        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomodel.getListOfSpecies()
        parameters_list = biomodel.getListOfParameters()
        reactions_list = biomodel.getListOfReactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = Species.getCurrentIndex()

        columns = Reaction.getCurrentIndex()

        self.forward_stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_reactants = individual_reaction.getListOfReactants()

            for individual_reactant in reaction_reactants:

                row = individual_reactant.index

                stoichiometry = individual_reactant.getStoichiometry()

                self.forward_stoichiometric_matrix[row, column] = int(stoichiometry)

        return self.forward_stoichiometric_matrix
    


    # ********************************
    # *           Function           *
    # ********************************
    def reverse_stoichiometric_matrix_constructor(self, biomodel) -> np.ndarray:
        """
        Constructs the reverse stoichiometric matrix for the given BioModel.

        Parameters:
            biomodel (an instance of BioModel class): The biological model containing species and reactions.

        Returns:
            np.ndarray: A 2D array representing the stoichiometric matrix, 
                        where rows correspond to species and columns to reactions.
        """

        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomodel.getListOfSpecies()
        parameters_list = biomodel.getListOfParameters()
        reactions_list = biomodel.getListOfReactions()

        if len(species_list) == 0:
            raise exceptions.EmptyList("There are no species in this model.")
        
        if len(reactions_list) == 0:
            raise exceptions.EmptyList("There are no reactions in this model.")

        rows = Species.getCurrentIndex()

        columns = Reaction.getCurrentIndex()

        self.reverse_stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            column = individual_reaction.index

            if column == None:
                continue

            reaction_products = individual_reaction.getListOfProducts()

            for individual_product in reaction_products:

                row = individual_product.index

                stoichiometry = individual_product.getStoichiometry()

                self.reverse_stoichiometric_matrix[row, column] = int(stoichiometry)

        return self.reverse_stoichiometric_matrix


    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_column_names(self, biomodel) -> dict:
        
        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        
        reactions_list = biomodel.getListOfReactions()

        column_indices_names = {}

        for reaction in reactions_list:

            column_indices_names[reaction.index] = reaction.ID

        return column_indices_names
    



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_row_names(self, biomodel) -> dict:
    
        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_list = biomodel.getListOfSpecies()

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
    def stoichiometric_matrix_element_information(self, i, j, biomodel, printing="off") -> str:

        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        highest_i = self.stoichiometric_matrix.shape[0]
        highest_j = self.stoichiometric_matrix.shape[1]

        if i > highest_i or j > highest_j:

            utility.error_printer("Invlaid indices provided!")
            
            print(f"\nThe highest value for i (rows), j (columns) are {highest_i} and {highest_j}, respectively")
            return

        else:

            row_indices_names = self.stoichiometric_matrix_row_names(biomodel)
            species = row_indices_names[i]

            column_indices_names = self.stoichiometric_matrix_column_names(biomodel)
            reaction = column_indices_names[j]

            if printing.lower() == "on":
                utility.printer(f"\nThe stoichiometric coefficient for {species} in reaction {reaction} is: ", f"{self.stoichiometric_matrix[i][j]}")

            return f"The stoichiometric coefficient for {species} in reaction {reaction} is: {self.stoichiometric_matrix[i][j]}"
    


    # ********************************
    # *           Function           *
    # ********************************
    def kinetic_constants_vector_constructor(self, biomodel, printing = "off") -> np.ndarray:

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        
        biomodel_reactions = biomodel.getListOfReactions()

        reactions_number = Reaction.getCurrentIndex()

        kinetic_vector_length = reactions_number

        vector_of_kinetic_constants = np.zeros(kinetic_vector_length)

        for biomodel_reaction in biomodel_reactions:

            if biomodel_reaction.index != None:

                index = biomodel_reaction.index

                name = biomodel_reaction.getId()

                k_plus_value = biomodel_reaction.kinetic_forward_rate_constant_value if biomodel_reaction.kinetic_forward_rate_constant_value is not None else 0.0

                k_minus_value = biomodel_reaction.kinetic_reverse_rate_constant_value if biomodel_reaction.kinetic_reverse_rate_constant_value is not None else 0.0

                if k_minus_value == 0.:
                    raise exceptions.NoReverseRateConstant(f"Kinetic Constants Vector cannot be constructed since there is no reverse reaction rate constant for reaction {name}")

                vector_of_kinetic_constants[index] = k_plus_value / k_minus_value

        if printing.lower() == "on":
            utility.printer("\nKinetic Constants Vector is:\n",vector_of_kinetic_constants)

        return vector_of_kinetic_constants
    

    # ********************************
    # *           Function           *
    # ********************************
    def kinetic_thermo_conversion_matrix_constructor(self, biomodel, printing = "off") -> np.ndarray:

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        reactions_number = Reaction.getCurrentIndex()

        identity_array = np.eye(reactions_number)

        forward_stoichiometric_matrix = self.forward_stoichiometric_matrix_constructor(biomodel)

        reverse_stoichiometric_matrix = self.reverse_stoichiometric_matrix_constructor(biomodel)

        transposed_forward_stoichiometric_matrix = np.transpose(forward_stoichiometric_matrix)

        transposed_reverse_stoichiometric_matrix = np.transpose(reverse_stoichiometric_matrix)

        conversion_matrix = np.block( [ [ identity_array, transposed_forward_stoichiometric_matrix ], [ identity_array, transposed_reverse_stoichiometric_matrix ] ] )

        if printing.lower() == "on":
            utility.printer("\nConversion Matrix is\n:", conversion_matrix)

        return conversion_matrix
    

    # ********************************
    # *           Function           *
    # ********************************
    def kinetic_rates_thermo_compatibility_check(self, biomodel, printing = "off") -> bool:

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        

        kinetic_rates_vector = self.kinetic_constants_vector_constructor(biomodel)

        with warnings.catch_warnings():

            warnings.simplefilter('error', RuntimeWarning)

            logn_kinetic_rates_vector = np.log(kinetic_rates_vector)

        logn_kinetic_rates_vector = logn_kinetic_rates_vector.reshape(-1,1)

        stoichiometric_matrix = self.stoichiometric_matrix_constructor(biomodel)

        minus_stoichiometric_matrix = -1 * stoichiometric_matrix

        minus_stoichiometric_null_space = null_space(minus_stoichiometric_matrix)

        transposed_minus_stoichiometric_null_space = minus_stoichiometric_null_space.T

        result = transposed_minus_stoichiometric_null_space @ logn_kinetic_rates_vector

        if np.all(np.abs(result) <= 1e-2):

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are compatible with thermodynamic constraints\n", text_color="green", text_style="bold")

            return True
        
        else:

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints\n", text_color="red", text_style="bold")
            
            return False