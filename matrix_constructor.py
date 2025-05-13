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

import itertools
    



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

            utility.error_printer("\nError: ", "Invlaid indices provided!")
            
            print(f"\nThe highest value for i (rows), j (columns) are {highest_i} and {highest_j}, respectively")
            return

        else:

            row_indices_names = self.stoichiometric_matrix_row_names(biomodel)
            species = row_indices_names[i]

            column_indices_names = self.stoichiometric_matrix_column_names(biomodel)
            reaction = column_indices_names[j]

            if printing.lower() == "on":
                utility.printer(f"\nThe stoichiometric coefficient for {species} in reaction {reaction} is: ", f"{self.stoichiometric_matrix[i][j]}", text_color="yellow")

            return f"The stoichiometric coefficient for {species} in reaction {reaction} is: {self.stoichiometric_matrix[i][j]}"


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
            rates = defaultdict(list)

            if expression.is_Mul:

                rates["forward_rate"].append(expression)

            elif expression.is_Add:

                for argument in expression.args:

                    if "-" in str(argument):

                        rates["reverse_rate"].append(argument)

                    else:

                        rates["forward_rate"].append(argument)

            return rates
        # --------------------------------------------------
        # --------------------------------------------------

        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def power_operator_finder(expression, parameters, parameters_values, new_parameters=None, reaction_rate_constant='', reaction_rate_constant_value=None, number_of_recursions = 0):

            # Create a dictionary to map operators to their corresponding operations
            operations = {
                "**": lambda x, y: pow(x, y)
            }

            if number_of_recursions > MAX_RECURSION:
                return expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value
            
            if new_parameters is None:
                new_parameters = []

            pairs = list(itertools.combinations(parameters, 2))

            used_elements = set()

            for pair in pairs:

                if pair[0] in used_elements or pair[1] in used_elements:
                    continue

                power_pattern = rf"(?:({re.escape(pair[0])})[\s]*(\*\*)[\s]*({re.escape(pair[1])})|({re.escape(pair[1])})[\s]*(\*\*)[\s]*({re.escape(pair[0])}))"

                power_match = re.search( power_pattern, expression )

                if power_match:

                    power_match_str = str(power_match.group(0))

                    used_elements.update(pair[0], pair[1])

                    parameters.remove(pair[0])

                    parameters.remove(pair[1])

                    replacing_power_k = "replaced_power_k_" + str(number_of_recursions)

                    parameters.append(replacing_power_k)

                    new_parameters.append(replacing_power_k)

                    if power_match.group(1):

                        if not bool(reaction_rate_constant):
                            reaction_rate_constant = str(power_match.group(1)) + str(power_match.group(2)) + str(power_match.group(3))
                        else:
                            reaction_rate_constant =  (reaction_rate_constant + str(power_match.group(2)) + str(power_match.group(3)) if str(power_match.group(1)) in new_parameters else str(power_match.group(1)) + str(power_match.group(2)) + reaction_rate_constant)

                        reaction_rate_constant_value = operations[str(power_match.group(2))](parameters_values[pair[0]], parameters_values[pair[1]])

                        parameters_values[replacing_power_k] = operations[str(power_match.group(2))](parameters_values[pair[0]], parameters_values[pair[1]])

                        expression = expression.replace( power_match_str, replacing_power_k )

                    elif power_match.group(4):

                        if not bool(reaction_rate_constant):
                            reaction_rate_constant = str(power_match.group(4)) + str(power_match.group(5)) + str(power_match.group(6))
                        else:
                            reaction_rate_constant =  (reaction_rate_constant + str(power_match.group(5)) + str(power_match.group(6)) if str(power_match.group(4)) in new_parameters else str(power_match.group(4)) + str(power_match.group(5)) + reaction_rate_constant)

                        reaction_rate_constant_value = operations[str(power_match.group(5))](parameters_values[pair[1]], parameters_values[pair[0]])

                        parameters_values[replacing_power_k] = operations[str(power_match.group(5))](parameters_values[pair[1]], parameters_values[pair[0]])

                        expression = expression.replace( power_match_str, replacing_power_k )

            return power_operator_finder(expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value, number_of_recursions+1)
        # --------------------------------------------------
        # --------------------------------------------------


        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def parameters_finder(expression, parameters, parameters_values, new_parameters = None, reaction_rate_constant = "", reaction_rate_constant_value=None, number_of_recursions=0):

            # Create a dictionary to map operators to their corresponding operations
            operations = {
                "*": lambda x, y: x * y,
                "+": lambda x, y: x + y,
                "-": lambda x, y: x - y,
                "/": lambda x, y: x / y,
            }

            if number_of_recursions > MAX_RECURSION:
                return expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value


            if new_parameters is None:
                new_parameters = []

            pairs = list(itertools.combinations(parameters, 2))

            used_elements = set()

            for pair in pairs:

                if pair[0] in used_elements or pair[1] in used_elements:
                    continue

                pattern = rf"(?:({re.escape(str(pair[0]))}).*?([\+\-\*/]).*?({re.escape(str(pair[1]))})|({re.escape(str(pair[1]))}).*?([\+\-\*/]).*?({re.escape(str(pair[0]))}))"

                match = re.search( pattern, str(expression) )

                if match:

                    match_str = str(match.group(0))

                    used_elements.update(pair[0], pair[1])

                    parameters.remove(pair[0])

                    parameters.remove(pair[1])

                    replacing_k = "replaced_k_" + str(number_of_recursions)

                    parameters.append(replacing_k)

                    new_parameters.append(replacing_k)

                    if match.group(1):

                        if not bool(reaction_rate_constant):
                            reaction_rate_constant = str(match.group(1)) + str(match.group(2)) + str(match.group(3))
                        else:
                            reaction_rate_constant =  (reaction_rate_constant + str(match.group(2)) + str(match.group(3)) if str(match.group(1)) in new_parameters else str(match.group(1)) + str(match.group(2)) + reaction_rate_constant)

                        reaction_rate_constant_value = operations[str(match.group(2))](parameters_values[pair[0]], parameters_values[pair[1]])

                        parameters_values[replacing_k] = operations[str(match.group(2))](parameters_values[pair[0]], parameters_values[pair[1]])

                        expression = expression.replace( match_str, replacing_k )

                    elif match.group(4):

                        if not bool(reaction_rate_constant):
                            reaction_rate_constant = str(match.group(4)) + str(match.group(5)) + str(match.group(6))
                        else:
                            reaction_rate_constant =  (reaction_rate_constant + str(match.group(5)) + str(match.group(6)) if str(match.group(4)) in new_parameters else str(match.group(4)) + str(match.group(5)) + reaction_rate_constant)

                        reaction_rate_constant_value = operations[str(match.group(5))](parameters_values[pair[1]], parameters_values[pair[0]])

                        parameters_values[replacing_k] = operations[str(match.group(5))](parameters_values[pair[1]], parameters_values[pair[0]])

                        expression = expression.replace( match_str, replacing_k )

            return parameters_finder(expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value, number_of_recursions+1)
        # --------------------------------------------------
        # --------------------------------------------------


        empty_global_parameters = False     # a flag to indicate if there are no global parameters

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")


        species_classes_list = biomodel.getListOfSpecies()
        parameter_classes_list = biomodel.getListOfParameters()
        reaction_classes_list = biomodel.getListOfReactions()
        function_definitions_list = biomodel.getListOfFunctionDefinitions()
        compartments_list = biomodel.getListOfCompartments()

        if not species_classes_list:
            raise ValueError("Species list is empty.")

        if not parameter_classes_list:
            empty_global_parameters = True

        if not reaction_classes_list:
            raise ValueError("Reaction list is empty.")



        string_to_sympy_symbols = {} # THERE ARE SOME BUILT-IN FUNCTIONS IN SYMPY THAT INTERFERE WITH VARIABLE NAMES AND THE STRING TO ... \
                                        # SYMPY EXPRESSIONS AND THE EXECUTION FAILS. TO PREVENT IT, I HAVE DEFINED SOME BUILT-IN FUNCTIONS AS ... \
                                        # SYMPY VARIABLES BEFOREHAND

        for individual_species_class in species_classes_list:
            species_name = individual_species_class.getId()
            string_to_sympy_symbols[species_name] = sp.symbols(species_name)

            string_to_sympy_symbols.update({
                name: sp.symbols(name) for name in compartments_list
            })


        parameters_values = {}  # This list stores the names of the parameters and their values

        for individual_parameter_class in parameter_classes_list:

            parameter_name = individual_parameter_class.getId()
            parameter_value = individual_parameter_class.getValue()
            string_to_sympy_symbols[parameter_name] = sp.symbols(parameter_name)
            parameters_values[parameter_name] = parameter_value


        for individual_reaction_class in reaction_classes_list:

            if individual_reaction_class.local_parameters is not None:

                local_parameters_values = {}

                local_parameters_list = individual_reaction_class.local_parameters

                for local_parameter_class in local_parameters_list:

                    local_parameter_name = local_parameter_class.getId()
                    local_parameter_value = local_parameter_class.getValue()
                    string_to_sympy_symbols[local_parameter_name] = sp.symbols(local_parameter_name)
                    local_parameters_values[local_parameter_name] = local_parameter_value

                parameters_values.update(local_parameters_values)
            
            else:

                if empty_global_parameters == True:
                    raise ValueError(f"There are no global parameters defined for the model, nor are there any local parameters defined for {reaction_name}!")

            reaction_name = individual_reaction_class.getId()
            reaction_rate_formula = individual_reaction_class.getKineticLaw()

            reaction_rate_formula, function_symbols = MatrixConstructor._expandFormula(reaction_rate_formula, function_definitions_list)

            string_to_sympy_symbols.update(function_symbols)

            try:
                sp_reaction_rate_formula = sp.sympify(reaction_rate_formula, locals = string_to_sympy_symbols, evaluate = False)

                simplified_formula = sp.cancel(sp_reaction_rate_formula)

                expanded_formula = sp.expand(simplified_formula)

            except sp.SympifyError as e:

                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")
                raise ValueError(f"\nSympify Error: {e}")

            except ValueError as v:
                utility.error_printer("\nValue Error: ", e)
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")
                raise ValueError(f"\nValue Error: {e}")

            except Exception as e:
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")
                raise ValueError(f"\nUnexpected Error: {e}")

            else:

                if printing.lower() == 'on':
                    utility.printer(f"\nThe simplified reaction rate expression for reaction {reaction_name} is:\n", simplified_formula, text_color="yellow")

            

            forward_reverse_rate_equations = get_forward_reverse_rate_expressions(expanded_formula)

            forward_rate_expressions = forward_reverse_rate_equations.get("forward_rate")

            if not forward_rate_expressions:
                raise ValueError(f"Forward rate expression cannot be found for {reaction_name}")
            
            
            forward_rate_expression = forward_rate_expressions[0]

            forward_variables_symbols = sp.sympify(forward_rate_expression, locals = string_to_sympy_symbols, evaluate = False).free_symbols
            forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
            forward_rate_matching_parameters = set(forward_variables_as_strings) & set(parameters_values.keys())

            
            if isinstance(forward_rate_matching_parameters, str):
                forward_rate_matching_parameters = [forward_rate_matching_parameters]
            else:
                forward_rate_matching_parameters = list(forward_rate_matching_parameters)

            forward_rate_constant = forward_rate_matching_parameters[0]
            forward_rate_constant_value = parameters_values[forward_rate_constant]

            if len( forward_rate_matching_parameters ) > 2 and len( forward_rate_matching_parameters ) < MAX_RECURSION:

                parameters_list = [str(p) for p in forward_rate_matching_parameters]

                expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value = power_operator_finder(str(forward_rate_expression), parameters_list, parameters_values)

                expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value = parameters_finder(expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value)

                if len(new_parameters) > 1:

                    message = f"\nThe forward kinetic rate constant for reaction {reaction_name} has more than one variable: {forward_rate_constant}"

                    utility.add_warning(message)
            
            elif len( forward_rate_matching_parameters ) >= MAX_RECURSION:

                raise ValueError(f"Number of reaction rate constants are more than {MAX_RECURSION} which is not supported!")
                

            for forward_rate_expression in forward_rate_expressions[1:]:

                forward_variables_symbols = sp.sympify(forward_rate_expression, locals = string_to_sympy_symbols, evaluate = False).free_symbols
                forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
                forward_rate_matching_parameters = set(forward_variables_as_strings) & set(parameters_values.keys())

                
                if isinstance(forward_rate_matching_parameters, str):
                    forward_rate_matching_parameters = [forward_rate_matching_parameters]
                else:
                    forward_rate_matching_parameters = list(forward_rate_matching_parameters)

                temp_forward_rate_constant = forward_rate_matching_parameters[0]
                temp_forward_rate_constant_value = parameters_values[temp_forward_rate_constant]

                if len( forward_rate_matching_parameters ) > 2 and len( forward_rate_matching_parameters ) <= MAX_RECURSION:

                    parameters_list = [str(p) for p in forward_rate_matching_parameters]

                    expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value = power_operator_finder(str(forward_rate_expression), parameters_list, parameters_values)

                    expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value = parameters_finder(expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value)

                    message = f"\nThe forward kinetic rate constant for reaction {reaction_name} has more than one variable: {forward_rate_constant}"

                    utility.add_warning(message)

                elif len( forward_rate_matching_parameters ) > MAX_RECURSION:

                    raise ValueError(f"Number of reaction rate constants are more than {MAX_RECURSION} which is not supported!")
                
                forward_rate_constant += " + " + temp_forward_rate_constant
                forward_rate_constant_value += temp_forward_rate_constant_value

                message = f"\nWARNING:\nThe forward kinetic rate constant for reaction {reaction_name} has more than one variable: {forward_rate_constant}"

                utility.message_printer(message, color="light_red", style ="normal")

                utility.add_warning(message)
            

            if forward_rate_constant:

                individual_reaction_class.kinetic_forward_rate_constant = forward_rate_constant
                individual_reaction_class.kinetic_forward_rate_constant_value = forward_rate_constant_value

                reverse_rate_constant = None

                if forward_reverse_rate_equations.get("reverse_rate"):

                    reverse_rate_expressions = forward_reverse_rate_equations.get("reverse_rate")

                    reverse_rate_expression = reverse_rate_expressions[0]

                    reverse_variables_symbols = sp.sympify(reverse_rate_expression, locals = string_to_sympy_symbols, evaluate = False).free_symbols
                
                    reverse_variables_as_strings = [str(symbol) for symbol in reverse_variables_symbols]

                    reverse_rate_matching_parameters = set(reverse_variables_as_strings) & set(parameters_values.keys())

                    
                    if isinstance(reverse_rate_matching_parameters, str):
                        reverse_rate_matching_parameters = [reverse_rate_matching_parameters]
                    else:
                        reverse_rate_matching_parameters = list(reverse_rate_matching_parameters)

                    reverse_rate_constant = reverse_rate_matching_parameters[0]
                    reverse_rate_constant_value = parameters_values[reverse_rate_constant]

                    if len( reverse_rate_matching_parameters ) > 2 and len( reverse_rate_matching_parameters ) <= MAX_RECURSION:

                        parameters_list = [str(p) for p in reverse_rate_matching_parameters]

                        expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value = power_operator_finder(str(reverse_rate_expression), parameters_list, parameters_values)

                        expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value = parameters_finder(expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value)

                        message = f"\nThe reverse kinetic rate constant for reaction {reaction_name} has more than one variable: {reverse_rate_constant}"

                        utility.add_warning(message)

                    elif len( reverse_rate_matching_parameters ) > MAX_RECURSION:

                        raise ValueError(f"Number of reaction rate constants are more than {MAX_RECURSION} which is not supported!")
                    


                    for reverse_rate_expression in reverse_rate_expressions[1:]:

                        reverse_variables_symbols = sp.sympify(reverse_rate_expression, locals = string_to_sympy_symbols, evaluate = False).free_symbols
                
                        reverse_variables_as_strings = [str(symbol) for symbol in reverse_variables_symbols]

                        reverse_rate_matching_parameters = set(reverse_variables_as_strings) & set(parameters_values.keys())

                        
                        if isinstance(reverse_rate_matching_parameters, str):
                            reverse_rate_matching_parameters = [reverse_rate_matching_parameters]
                        else:
                            reverse_rate_matching_parameters = list(reverse_rate_matching_parameters)

                        temp_reverse_rate_constant = reverse_rate_matching_parameters[0]
                        temp_reverse_rate_constant_value = parameters_values[temp_reverse_rate_constant]
                        
                        
                        if len( reverse_rate_matching_parameters ) > 2 and len( reverse_rate_matching_parameters ) <= MAX_RECURSION:

                            parameters_list = [str(p) for p in reverse_rate_matching_parameters]

                            expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value = power_operator_finder(str(reverse_rate_expression), parameters_list, parameters_values)

                            expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value = parameters_finder(expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value)

                            message = f"\nThe reverse kinetic rate constant for reaction {reaction_name} has more than one variable: {reverse_rate_constant}"

                            utility.add_warning(message)
                        
                        elif len( reverse_rate_matching_parameters ) > MAX_RECURSION:

                            raise ValueError(f"Number of reaction rate constants are more than {MAX_RECURSION} which is not supported!")
                                
                        reverse_rate_constant += " + " + temp_reverse_rate_constant
                        reverse_rate_constant_value += temp_reverse_rate_constant_value

                        message = f"\nThe reverse kinetic rate constant for reaction {reaction_name} has more than one variable: {reverse_rate_constant}"

                        utility.add_warning(message)
                        

                    individual_reaction_class.kinetic_reverse_rate_constant = reverse_rate_constant
                    individual_reaction_class.kinetic_reverse_rate_constant_value = reverse_rate_constant_value



                if forward_rate_constant is not None:

                    if reverse_rate_constant is None:

                        if printing.lower() == "on":

                            utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")

                        reaction_to_rate_constants[reaction_name] = {"forward": forward_rate_constant_value}

                    else:

                        if printing.lower() == "on":

                            utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")
                            utility.printer("Reverse rate constant is:", reverse_rate_constant, text_style="bold")


                        reaction_to_rate_constants[reaction_name] = {"forward": forward_rate_constant_value, "reverse": reverse_rate_constant_value}

            else:

                if forward_reverse_rate_equations.get("reverse_rate"):

                    forward_variables_symbols = sp.sympify(simplified_formula, locals = string_to_sympy_symbols, evaluate = False).free_symbols
                    forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
                    common_rate_constant = next(iter(set(forward_variables_as_strings) & set(parameters_values.keys())), None)

                    forward_rate_constant = common_rate_constant

                    reverse_rate_constant = common_rate_constant


                    individual_reaction_class.kinetic_forward_rate_constant = forward_rate_constant
                    individual_reaction_class.kinetic_forward_rate_constant_value = parameters_values[str(forward_rate_constant)]

                    individual_reaction_class.kinetic_reverse_rate_constant = reverse_rate_constant
                    individual_reaction_class.kinetic_reverse_rate_constant_value = parameters_values[str(reverse_rate_constant)]

                    if forward_rate_constant is not None:

                        if reverse_rate_constant is None:

                            if printing.lower() == "on":

                                utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")

                            reaction_to_rate_constants[reaction_name] = {"forward": parameters_values[str(forward_rate_constant)]}

                        else:

                            if printing.lower() == "on":

                                utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")
                                utility.printer("Reverse rate constant is:", reverse_rate_constant, text_style="bold")


                            reaction_to_rate_constants[reaction_name] = {"forward": parameters_values[str(forward_rate_constant)], "reverse": parameters_values[str(reverse_rate_constant)]}



        return reaction_to_rate_constants
    


    # ********************************
    # *           Function           *
    # ********************************
    def kinetic_constants_vector_constructor(self, biomodel, printing = "off") -> np.ndarray:

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        
        biomodel_reactions = biomodel.getListOfReactions()

        called = False

        for _ in range(3):

            r_num = random.randint(0, Reaction.getCurrentIndex()-1)

            if biomodel_reactions[r_num].kinetic_forward_rate_constant_value is not None:

                called = True

        if called == False:
        
            self.forward_reverse_rate_finder(biomodel, printing="on")

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
            utility.printer("\nKinetic Constants Vector is:\n",vector_of_kinetic_constants, text_color="green")

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
            utility.printer("\nConversion Matrix is\n:",conversion_matrix, text_color="magenta", text_style="bold")

        return conversion_matrix
    

    # ********************************
    # *           Function           *
    # ********************************
    def kinetic_rates_thermo_compatibility_check(self, biomodel, printing = "off") -> bool:

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")
        
        try:

            kinetic_rates_vector = self.kinetic_constants_vector_constructor(biomodel)

        except exceptions.NoReverseRateConstant as e:

            utility.error_printer(f"\nCaught an error: ", e)

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints", text_color="red", text_style="bold")

                return False

        with warnings.catch_warnings():

            warnings.simplefilter('error', RuntimeWarning)

            try:

                logn_kinetic_rates_vector = np.log(kinetic_rates_vector)

            except ValueError as e:
                utility.printer(f"\nCaught an error: ", e)
                print(f"\nCaught an error: {e}")

                if printing.lower() == "on":
                    utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints", text_color="red", text_style="bold")

                raise ValueError(f"\nCaught an error: ", e)

            except RuntimeWarning as e:

                utility.error_printer(f"\nCaught an error: ", e)

                if printing.lower() == "on":
                    utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints", text_color="red", text_style="bold")

                return False

        logn_kinetic_rates_vector = logn_kinetic_rates_vector.reshape(-1,1)

        stoichiometric_matrix = self.stoichiometric_matrix_constructor(biomodel)

        minus_stoichiometric_matrix = -1 * stoichiometric_matrix

        minus_stoichiometric_null_space = null_space(minus_stoichiometric_matrix)

        transposed_minus_stoichiometric_null_space = minus_stoichiometric_null_space.T

        result = transposed_minus_stoichiometric_null_space @ logn_kinetic_rates_vector

        if np.all(np.abs(result) <= 1e-2):

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are compatible with thermodynamic constraints", text_color="green", text_style="bold")

            return True
        
        else:

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints", text_color="red", text_style="bold")
            
            return False
        



    @staticmethod
    def _expandFormula(formula, function_definitions,
            num_recursions=0):
        """
        Expands the kinetics formula, replacing function definitions
        with their body.

        Parameters
        ----------
        expansion: str
            expansion of the kinetic law
        function_definitions: list-FunctionDefinition
        num_recursion: int
        
        Returns
        -------
        str
        """


        if num_recursions > MAX_RECURSION:
            return formula, {}
        
        input_symbols_dict = {}
        done = True
        for function_definition in function_definitions:
        
            function_calls = re.findall(r'{}\(.*?\)'.format(function_definition.ID), formula)

            if len(function_calls) == 0:
                continue

            done = False
            for function_call in function_calls:
                
                call_arguments = re.findall(r'\(.*?\)', function_call)[0]
                call_arguments = call_arguments.strip()
                call_arguments = call_arguments[1:-1]
                input_arguments = call_arguments.split(',')
                input_arguments = [a.strip() for a in input_arguments]
                for arg in input_arguments:
                    input_symbols_dict[arg] = sp.symbols(arg)
                function_formula = function_definition.formula
                for formal_arg, call_arg in zip(function_definition.arguments, input_arguments):
                    function_formula = function_formula.replace(formal_arg, call_arg)
                formula = formula.replace(function_call, function_formula)
        if not done:
            return MatrixConstructor._expandFormula(formula, function_definitions,
                        num_recursions=num_recursions+1)
        return formula, input_symbols_dict