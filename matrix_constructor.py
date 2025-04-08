import numpy as np
import exceptions
import os
import sympy as sp
from classes.cReaction import *
from classes.cSpecies import *
import utility

from colorama import Fore, Back, Style, init

    



class MatrixConstructor:



    # # ********************************
    # # *           Function           *
    # # ********************************
    # def stoichiometric_matrix_constructor(self, biomodel) -> np.ndarray:
    #     """
    #     Constructs the stoichiometric matrix for the given BioModel.

    #     Parameters:
    #         biomodel (an instance of BioModel class): The biological model containing species and reactions.

    #     Returns:
    #         np.ndarray: A 2D array representing the stoichiometric matrix, 
    #                     where rows correspond to species and columns to reactions.
    #     """

    #     if biomodel == None:
    #         raise exceptions.NoModel("There is no input model")

    #     species_list = biomodel.getListOfSpecies()
    #     parameters_list = biomodel.getListOfParameters()
    #     reactions_list = biomodel.getListOfReactions()

    #     try:

    #         if len(species_list) == 0:
    #             raise exceptions.EmptyList("There are no species in this model.")
            
    #         if len(reactions_list) == 0:
    #             raise exceptions.EmptyList("There are no reactions in this model.")
            
    #     except exceptions.EmptyList as e:
    #         print(Fore.RED +  f"\nError: {e}")
    #         print("Unable to complete the query\!")
    #         return


    #     self.species_indices = {}

    #     current_element_index = 0

    #     for individual_species in species_list:

    #         species_name = individual_species.getId() # I AM NOT SURE IF I NEED TO USE getName() or getID() to READ THE NAMES OF SPECIES

    #         if species_name != "empty":
    #             if species_name not in self.species_indices:
    #                 self.species_indices[species_name] = current_element_index
    #                 current_element_index += 1


    #     self.reaction_indices = {}

    #     current_reaction_index = 0

    #     rows = len(self.species_indices)

    #     columns = len(reactions_list)

    #     self.stoichiometric_matrix = np.zeros((rows, columns), dtype = int)

    #     for individual_reaction in reactions_list:

    #         reaction_name = individual_reaction.getId()

    #         self.reaction_indices[reaction_name] = current_reaction_index
    #         column = current_reaction_index
    #         current_reaction_index += 1

    #         reaction_reactants = individual_reaction.getListOfReactants()

    #         for individual_reactant in reaction_reactants:
    #             reactant_name = individual_reactant.getId()
    #             reactant_stoichiometry = individual_reactant.getStoichiometry()

    #             if reactant_name != "empty":

    #                 row = self.species_indices[reactant_name]
    #                 self.stoichiometric_matrix[row, column] = -1 * int(reactant_stoichiometry)

    #         reaction_products = individual_reaction.getListOfProducts()

    #         if reactant_name == "empty":

    #             for individual_product in reaction_products:
    #                 product_name = individual_product.getSpecies()
    #                 product_stoichiometry = individual_product.getStoichiometry()

    #             row = self.species_indices[product_name]
    #             self.stoichiometric_matrix[row, column] = int(product_stoichiometry)

    #             new_reactant = product_name + "_e"
    #             self.species_indices[new_reactant] = current_element_index
    #             row = current_element_index
    #             current_element_index += 1

    #             new_row = np.zeros((1, self.stoichiometric_matrix.shape[1]))
    #             self.stoichiometric_matrix = np.vstack([self.stoichiometric_matrix, new_row])

    #             self.stoichiometric_matrix[row, column] = -1

    #         else:

    #             for individual_product in reaction_products:
    #                 product_name = individual_product.getSpecies()
    #                 product_stoichiometry = individual_product.getStoichiometry()

    #                 if product_name == "empty":
    #                     new_product = reactant_name + "_e"
    #                     self.species_indices[new_product] = current_element_index
    #                     row = current_element_index
    #                     current_element_index += 1

    #                     new_row = np.zeros((1, self.stoichiometric_matrix.shape[1]))
    #                     self.stoichiometric_matrix = np.vstack([self.stoichiometric_matrix, new_row])

    #                     self.stoichiometric_matrix[row, column] = 1

    #                 else:

    #                     row = self.species_indices[product_name]
    #                     self.stoichiometric_matrix[row, column] = int( product_stoichiometry )

    #     return self.stoichiometric_matrix


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
            raise exceptions.NoModel("There is no input model")

        species_list = biomodel.getListOfSpecies()
        parameters_list = biomodel.getListOfParameters()
        reactions_list = biomodel.getListOfReactions()

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            print(Fore.RED +  f"\nError: {e}")
            print("Unable to complete the query\!")
            return

        rows = Species.getCurrentIndex() + 1

        columns = Reaction.getCurrentIndex() + 1

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
            raise exceptions.NoModel("There is no input model")

        species_list = biomodel.getListOfSpecies()
        parameters_list = biomodel.getListOfParameters()
        reactions_list = biomodel.getListOfReactions()

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            print(Fore.RED +  f"\nError: {e}")
            print("Unable to complete the query\!")
            return

        rows = Species.getCurrentIndex() + 1

        columns = Reaction.getCurrentIndex() + 1

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
            raise exceptions.NoModel("There is no input model")

        species_list = biomodel.getListOfSpecies()
        parameters_list = biomodel.getListOfParameters()
        reactions_list = biomodel.getListOfReactions()

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            print(Fore.RED +  f"\nError: {e}")
            print("Unable to complete the query\!")
            return

        rows = Species.getCurrentIndex() + 1

        columns = Reaction.getCurrentIndex() + 1

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
    def stoichiometric_matrix_column_names(self) -> dict:
        
        try:
            self.reaction_indices

            return self.reaction_indices
        
        except NameError:
            print(Fore.RED + "Stoichiometric Matrix hasn't been instantiated yet!\nPlease call the Stoichiometric Matrix Constructor Function first.")
    



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_row_names(self) -> dict:
    
        try:
            self.species_indices

            return self.species_indices
        
        except NameError:
            print(Fore.RED + "Stoichiometric Matrix hasn't been instantiated yet!\nPlease call the Stoichiometric Matrix Constructor Function first.")
    



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_element_information(self, i, j) -> None:

        reaction = next((k for k, v in self.reaction_indices.items() if v == j), None)

        species = next((k for k, v in self.species_indices.items() if v == i), None)

        if (reaction is None) or (species is None):
            print(Fore.RED + "Error: Invlaid indices provided!")
            highest_i = self.stoichiometric_matrix.shape[0]
            highest_j = self.stoichiometric_matrix.shape[1]
            print(Fore.MAGENTA + f"\nThe highest value for i (rows), j (columns) are {highest_i} and {highest_j}, respectively")
            return

        print(Fore.YELLOW + f"The stoichiometric coefficient for {species} in reaction {reaction} is {self.stoichiometric_matrix[i][j]}")


    # ********************************
    # *           Function           *
    # ********************************
    def forward_reverse_rate_finder(self, biomodel, printing = "off") -> dict:
        """
        Input of this function is a BioModel class instance
        This function extracts the rate constants for each reaction in the given model, 
        including both forward and reverse directions. Then, "kinetic_forward_rate_constant" and "kinetic_reverse_rate_constant" variables of each reaction is updated

        Each reaction's rate constants are stored in a dictionary with 'forward' and 'reverse' 
        as keys and their corresponding values as floats.

        These per-reaction dictionaries are then stored in a master dictionary, 
        `reaction_to_rate_constants`, where the keys are reaction names and the values 
        are the corresponding rate constant dictionaries.

        The final dictionary is returned.
        """

        reaction_to_rate_constants = {} # This dictionary maps reaction names (IDs) to a dictionary of forward and reverse rate constants and their values

        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def get_forward_reverse_rate_expressions(expression):

            # Separate positive and negative terms manually
            rates = {}

            if "-" in str(expression):

                if expression.is_Mul:

                    for arg in expression.args:

                        if arg.is_Add:

                            for term in arg.args:
                    
                                if "-" in str(term):
                                    rates["reverse_rate"] = term
                                    
                                else:
                                    rates["forward_rate"] = term

                        elif arg.is_Mul:
                            nested_rates = get_forward_reverse_rate_expressions(arg)
                            rates.update(nested_rates)

                elif expression.is_Add:

                    for arg in expression.args:

                        if arg.is_Add:

                            for term in arg.args:
                    
                                if "-" in str(term):
                                    rates["reverse_rate"] = term
                                    
                                else:
                                    rates["forward_rate"] = term

            else:
                rates["forward_rate"] = expression

            return rates
        # --------------------------------------------------
        # --------------------------------------------------


        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_classes_list = biomodel.getListOfSpecies()
        parameter_classes_list = biomodel.getListOfParameters()
        reaction_classes_list = biomodel.getListOfReactions()

        string_to_sympy_symbols = {} # THERE ARE SOME BUILT-IN FUNCTIONS IN SYMPY THAT INTERFERES WITH VARIABLE NAMES AND THE STRING TO ... \
                                        # SYMPY EXPRESSIONS AND THE EXECUTION FAILS. TO PREVENT IT, I HAVE DEFINED SOME BUILT-IN FUNCTIONS AS ... \
                                        # SYMPY VARIABLES BEFOREHAND

        for individual_species_class in species_classes_list:
            species_name = individual_species_class.getId()
            string_to_sympy_symbols[species_name] = sp.symbols(species_name)


        parameters_values = {}  # This list stores the names of the parameters and their values

        for individual_parameter_class in parameter_classes_list:

            parameter_name = individual_parameter_class.getId()
            parameter_value = individual_parameter_class.getValue()
            string_to_sympy_symbols[parameter_name] = sp.symbols(parameter_name)
            parameters_values[parameter_name] = parameter_value


        for individual_reaction_class in reaction_classes_list:

            reaction_name = individual_reaction_class.getId()
            reaction_rate_formula = individual_reaction_class.getKineticLaw()


            try:
                sp_reaction_rate_formula = sp.sympify(reaction_rate_formula, locals = string_to_sympy_symbols, evaluate = False)

            except sp.SympifyError as e:
                print(Fore.RED + f"\nSympify Error: {e}" + Fore.LIGHTRED_EX + "Function cannot be completed")
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")

            except ValueError as v:
                print(f"\nValue Error: {v}" + Fore.LIGHTRED_EX + "Function cannot be completed")
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")

            except Exception as e:
                print(f"\nUnexpected Error: {e}" + Fore.LIGHTRED_EX + "Function cannot be completed")
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")

            else:

                if printing.lower() == 'on':
                    utility.printer("\nThe reaction rate expression is:\n", sp_reaction_rate_formula, text_color="yellow")
                #print(f"\nThe reaction rate expression is:" + Fore.YELLOW + f"\n{sp_reaction_rate_formula}")

            simplified_formula = sp.simplify(sp_reaction_rate_formula)

            forward_reverse_rate_equations = get_forward_reverse_rate_expressions(simplified_formula)

            forward_variables_symbols = sp.sympify(forward_reverse_rate_equations.get("forward_rate"), locals = string_to_sympy_symbols, evaluate = False).free_symbols
            forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
            forward_rate_constant = next(iter(set(forward_variables_as_strings) & set(parameters_values.keys())), None)

            individual_reaction_class.kinetic_forward_rate_constant = forward_rate_constant
            individual_reaction_class.kinetic_forward_rate_constant_value = parameters_values[str(forward_rate_constant)]

            reverse_rate_constant = None

            if forward_reverse_rate_equations.get("reverse_rate"):

                reverse_variables_symbols = sp.sympify(forward_reverse_rate_equations.get("reverse_rate"), locals = string_to_sympy_symbols, evaluate = False).free_symbols
            
                reverse_variables_as_strings = [str(symbol) for symbol in reverse_variables_symbols]

                reverse_rate_constant = next(iter(set(reverse_variables_as_strings) & set(parameters_values.keys())), None)

                individual_reaction_class.kinetic_reverse_rate_constant = reverse_rate_constant
                individual_reaction_class.kinetic_reverse_rate_constant_value = parameters_values[str(reverse_rate_constant)]


            if reverse_rate_constant is None:

                if printing.lower() == "on":

                    utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")

                #print(f"\nForward rate constant is " + Style.BRIGHT + Fore.CYAN + f"{forward_rate_constant}")

                reaction_to_rate_constants[reaction_name] = {"forward": parameters_values[str(forward_rate_constant)]}

            else:

                if printing == "On" or printing == "ON" or printing == "on":

                    utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")
                    utility.printer("Reverse rate constant is:", reverse_rate_constant, text_style="bold")


                reaction_to_rate_constants[reaction_name] = {"forward": parameters_values[str(forward_rate_constant)], "reverse": parameters_values[str(reverse_rate_constant)]}

        return reaction_to_rate_constants
    


    # ********************************
    # *           Function           *
    # ********************************
    # def conversion_matrix_constructor(self, biomodel, printing = "off"):

