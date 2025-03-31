import numpy as np
import exception
import os




class MatrixConstructor:

    def SBML_stoichiomrtic_matrix_constructor(self, model):

        species_list = model.getListOfSpecies()
        parameters_list = model.getListOfParameters()
        reactions_list = model.getListOfReactions()
        rules_list = model.getListOfRules() #The governing rules of the reactions

        try:

            if len(species_list) == 0:
                raise exception.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exception.EmptyList("There are no reactions in this model.")
            
        except exception.EmptyList as e:
            print(f"\nError: {e}")


        self.species_indices = {}

        current_element_index = 0

        for individual_species in species_list:

            species_name = individual_species.getId() # I AM NOT SURE IF I NEED TO USE getName() or getID() to READ THE NAMES OF SPECIES

            if species_name != "empty":
                if species_name not in self.species_indices:
                    self.species_indices[species_name] = current_element_index
                    current_element_index += 1

        self.reaction_indices = {}

        current_reaction_index = 0

        rows = len(self.species_indices)

        columns = len(reactions_list)

        self.stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

        for individual_reaction in reactions_list:

            reaction_name = individual_reaction.getId()

            self.reaction_indices[reaction_name] = current_reaction_index
            column = current_reaction_index
            current_reaction_index += 1

            reaction_reactants = individual_reaction.getListOfReactants()

            for individual_reactant in reaction_reactants:
                reactant_name = individual_reactant.getSpecies()
                reactant_stoichiometry = individual_reactant.getStoichiometry()

                if reactant_name != "empty":

                    row = self.species_indices[reactant_name]
                    self.stoichiometric_matrix[row, column] = -1 * int(reactant_stoichiometry)

            reaction_products = individual_reaction.getListOfProducts()

            if reactant_name == "empty":

                for individual_product in reaction_products:
                    product_name = individual_product.getSpecies()
                    product_stoichiometry = individual_product.getStoichiometry()

                row = self.species_indices[product_name]
                self.stoichiometric_matrix[row, column] = int(product_stoichiometry)

                new_reactant = product_name + "_e"
                self.species_indices[new_reactant] = current_element_index
                row = current_element_index
                current_element_index += 1

                new_row = np.zeros((1, self.stoichiometric_matrix.shape[1]))
                self.stoichiometric_matrix = np.vstack([self.stoichiometric_matrix, new_row])

                self.stoichiometric_matrix[row, column] = -1

            else:

                for individual_product in reaction_products:
                    product_name = individual_product.getSpecies()
                    product_stoichiometry = individual_product.getStoichiometry()

                    if product_name == "empty":
                        new_product = reactant_name + "_e"
                        self.species_indices[new_product] = current_element_index
                        row = current_element_index
                        current_element_index += 1

                        new_row = np.zeros((1, self.stoichiometric_matrix.shape[1]))
                        self.stoichiometric_matrix = np.vstack([self.stoichiometric_matrix, new_row])

                        self.stoichiometric_matrix[row, column] = 1

                    else:

                        row = self.species_indices[product_name]
                        self.stoichiometric_matrix[row, column] = int( product_stoichiometry )

        return self.stoichiometric_matrix
    
    def stoichiometric_matrix_column_names(self):
        
        return self.reaction_indices
    
    def stoichiometric_matrix_row_names(self):

        return self.species_indices
    
    def elementinformation(self, i, j):

        reaction = next((k for k, v in self.reaction_indices.items() if v == j), None)

        species = next((k for k, v in self.species_indices.items() if v == i), None)

        if (reaction is None) or (species is None):
            print("Error: Invlaid indoces provided!")
            highest_i = self.stoichiometric_matrix.shape[0]
            highest_j = self.stoichiometric_matrix.shape[1]
            print(f"\nThe highest value for i (rows), j (columns) are {highest_i} and {highest_j}, respectively")
            return

        info = f"The stoichiometric coefficient for {species} in reaction {reaction} is {self.stoichiometric_matrix[i][j]}"

        print(f"The stoichiometric coefficient for {species} in reaction {reaction} is {self.stoichiometric_matrix[i][j]}")