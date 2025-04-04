import numpy as np
import exceptions
import os
import sympy as sp
from classes.cReaction import *

from colorama import Fore, Back, Style, init

    



class MatrixConstructor:



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiomrtic_matrix_constructor(self, model):
        '''
        This function constructs the Stoichiometric Matrix for a SBML model
        It takes a SBML model, which is already read by another function, and extracts species and reactions to construct the Stoichiometric Matrix
        '''

        if model == None:
            raise exceptions.NoModel("There is no input model")

        species_list = model.getListOfSpecies()
        parameters_list = model.getListOfParameters()
        reactions_list = model.getListOfReactions()
        # rules_list = model.getListOfRules() #The governing rules of the reactions

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            print(Fore.RED +  f"\nError: {e}")


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
                reactant_name = individual_reactant.getId()
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
    



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_column_names(self):
        
        try:
            self.reaction_indices

            return self.reaction_indices
        
        except NameError:
            print(Fore.RED + "Stoichiometric Matrix hasn't been instantiated yet!\nPlease call the Stoichiometric Matrix Constructor Function first.")
    



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_row_names(self):
    
        try:
            self.species_indices

            return self.species_indices
        
        except NameError:
            print(Fore.RED + "Stoichiometric Matrix hasn't been instantiated yet!\nPlease call the Stoichiometric Matrix Constructor Function first.")
    



    # ********************************
    # *           Function           *
    # ********************************
    def stoichiometric_matrix_element_information(self, i, j):

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
    def conversion_matrix_constructor(self, model):
        '''
        This function finds the reaction rate constants that are needed to check thermodynamic consistency
        '''

        # ################################
        # *      Internal Function       *
        # ################################
        def get_forward_reverse_rates(expression):

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
                            get_forward_reverse_rates(arg)

            else:
                rates["forward_rate"] = expression

            return rates





        init( autoreset=True )

        if model is None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        species_classes_list = model.getListOfSpecies()
        parameter_classes_list = model.getListOfParameters()
        reaction_classes_list = model.getListOfReactions()
        #rule_classes_list = model.getListOfRules() #The governing rules of the reactions
        parameters_list = []

        for individual_parameter_class in parameter_classes_list:

            parameter_name = individual_parameter_class.getId()
            parameters_list.append(parameter_name)
            # parameter_SBO_term_URL = individual_parameter_class.getSBOTermAsURL()
            # parameter_SBO_term = os.path.basename(parameter_SBO_term_URL)

        for individual_reaction_class in reaction_classes_list:

            reaction_name = individual_reaction_class.getId()
            reaction_rate_formula = individual_reaction_class.getKineticLaw()
            # kinetic_law_class = individual_reaction_class.getKineticLaw()
            # reaction_rate_formula = kinetic_law_class.getFormula()

            symbols_dict = {"alpha": sp.symbols("alpha"), "beta": sp.symbols("beta")} # THERE ARE SOME BUILT-IN FUNCTIONS IN SYMPY THAT INTERFERES WITH VARIABLE NAMES AND THE STRING TO ... \
                                                                                    # SYMPY EXPRESSIONS AND THE EXECUTION FAILS. TO PREVENT IT, I HAVE DEFINED SOME BUILT-IN FUNCTIONS AS ... \
                                                                                    # SYMPY VARIABLES BEFOREHAND

            try:
                sp_reaction_rate_formula = sp.sympify(reaction_rate_formula, locals = symbols_dict, evaluate = False)

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
                print(f"\nThe reaction rate expression is:" + Fore.YELLOW + f"\n{sp_reaction_rate_formula}")

            simplified_formula = sp.simplify(sp_reaction_rate_formula)

            reordered_simplified_formulae = get_forward_reverse_rates(simplified_formula)

            forward_variables_symbols = sp.sympify(reordered_simplified_formulae.get("forward_rate"), locals = symbols_dict, evaluate = False).free_symbols
            forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
            forward_rate_constat = set(forward_variables_as_strings) & set(parameters_list)

            reverse_rate_constant = None

            if reordered_simplified_formulae.get("reverse_rate"):

                reverse_variables_symbols = sp.sympify(reordered_simplified_formulae.get("reverse_rate"), locals = symbols_dict, evaluate = False).free_symbols
            
                reverse_variables_as_strings = [str(symbol) for symbol in reverse_variables_symbols]

                reverse_rate_constant = set(reverse_variables_as_strings) & set(parameters_list)


            if reverse_rate_constant is None:

                print(f"\nForward rate constant is " + Style.BRIGHT + Fore.CYAN + f"{list(forward_rate_constat)[0]}")

            else:

                print(f"\nForward rate constant is " + Style.BRIGHT + Fore.CYAN + f"{list(forward_rate_constat)[0]}" + Style.NORMAL + Fore.WHITE + " and reverse rate constant is " + Style.BRIGHT + Fore.CYAN + f"{list(reverse_rate_constant)[0]}")