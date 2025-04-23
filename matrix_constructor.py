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

            if reaction.index:

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

            if species.index:

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

                        elif arg.is_Mul:

                            if "-" in str(arg):
                                rates["reverse_rate"] = arg
                                    
                            else:
                                rates["forward_rate"] = arg

            else:
                rates["forward_rate"] = expression

            return rates
        # --------------------------------------------------
        # --------------------------------------------------


        empty_global_parameters = False     # a flag to indicate if there are no global parameters

        if biomodel is None:
            raise exceptions.NoModel("No BioModel has been read!!!")


        species_classes_list = biomodel.getListOfSpecies()
        parameter_classes_list = biomodel.getListOfParameters()
        reaction_classes_list = biomodel.getListOfReactions()

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


        parameters_values = {}  # This list stores the names of the parameters and their values

        for individual_parameter_class in parameter_classes_list:

            parameter_name = individual_parameter_class.getId()
            parameter_value = individual_parameter_class.getValue()
            string_to_sympy_symbols[parameter_name] = sp.symbols(parameter_name)
            parameters_values[parameter_name] = parameter_value


        for individual_reaction_class in reaction_classes_list:

            reaction_name = individual_reaction_class.getId()
            reaction_rate_formula = individual_reaction_class.getKineticLaw()

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

            try:
                sp_reaction_rate_formula = sp.sympify(reaction_rate_formula, locals = string_to_sympy_symbols, evaluate = False)

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
                    utility.printer("\nThe reaction rate expression is:\n", sp_reaction_rate_formula, text_color="yellow")

            simplified_formula = sp.simplify(sp_reaction_rate_formula)

            forward_reverse_rate_equations = get_forward_reverse_rate_expressions(simplified_formula)

            forward_rate_expr = forward_reverse_rate_equations.get("forward_rate")

            if forward_rate_expr is None:
                raise ValueError(f"Forward rate expression cannot be found for {reaction_name}")
            
            forward_variables_symbols = sp.sympify(forward_rate_expr, locals = string_to_sympy_symbols, evaluate = False).free_symbols
            forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
            forward_rate_matching_parameters = set(forward_variables_as_strings) & set(parameters_values.keys())

            if not forward_rate_matching_parameters:
                raise exceptions.EmptyList(f"No forward reaction rate constant found in reaction {reaction_name}")
            else:
                if isinstance(forward_rate_matching_parameters, str):
                    forward_rate_matching_parameters = [forward_rate_matching_parameters]
                else:
                    forward_rate_matching_parameters = list(forward_rate_matching_parameters)

                forward_rate_constant = forward_rate_matching_parameters[0]
                forward_rate_constant_value = parameters_values[forward_rate_constant]

                if len( forward_rate_matching_parameters ) == 2:

                    k_plus_1 = str(forward_rate_matching_parameters[0])

                    k_plus_2 = str(forward_rate_matching_parameters[1])

                    forward_pattern = rf"(?:{re.escape(k_plus_1)}[\w\s]*([\+\-\*/])[\w\s]*{re.escape(k_plus_2)}|{re.escape(k_plus_2)}[\w\s]*([\+\-\*/])[\w\s]*{re.escape(k_plus_1)})"

                    forward_match = re.search( forward_pattern, str(forward_rate_expr) )

                    forward_match_str = forward_match.group()

                    if "*" in forward_match_str:
                        forward_rate_constant = k_plus_1 + " * " + k_plus_2
                        forward_rate_constant_value = parameters_values[k_plus_1] *  parameters_values[k_plus_2]

                    elif "+" in forward_match_str:
                        forward_rate_constant = k_plus_1 + " + " + k_plus_2
                        forward_rate_constant_value = parameters_values[k_plus_1] +  parameters_values[k_plus_2]

                    elif "-" in forward_match_str:
                        forward_rate_constant = k_plus_1 + " - " + k_plus_2
                        forward_rate_constant_value = parameters_values[k_plus_1] -  parameters_values[k_plus_2]

                    elif "/" in forward_match_str:
                        forward_rate_constant = k_plus_1 + " / " + k_plus_2
                        forward_rate_constant_value = parameters_values[k_plus_1] /  parameters_values[k_plus_2]

                    else:
                        
                        raise ValueError(f"Operand not found in the reaction rate composition expression: {forward_rate_expr} in line {inspect.currentframe().f_lineno()} in \"matrix_constructor.py\"")
                    
                elif len( forward_rate_matching_parameters ) > 2:

                    raise ValueError("Number of reaction rate constants are more than two which is not supported!")
            

            if forward_rate_constant:

                individual_reaction_class.kinetic_forward_rate_constant = forward_rate_constant
                individual_reaction_class.kinetic_forward_rate_constant_value = forward_rate_constant_value

                reverse_rate_constant = None

                if forward_reverse_rate_equations.get("reverse_rate"):

                    reverse_rate_expr = forward_reverse_rate_equations.get("reverse_rate")

                    reverse_variables_symbols = sp.sympify(reverse_rate_expr, locals = string_to_sympy_symbols, evaluate = False).free_symbols
                
                    reverse_variables_as_strings = [str(symbol) for symbol in reverse_variables_symbols]

                    reverse_rate_matching_parameters = next(iter(set(reverse_variables_as_strings) & set(parameters_values.keys())), None)

                    if not reverse_rate_matching_parameters:
                        raise exceptions.EmptyList(f"No reverse reaction rate constant found in reaction {reaction_name}")
                    else:
                        if isinstance(reverse_rate_matching_parameters, str):
                            reverse_rate_matching_parameters = [reverse_rate_matching_parameters]
                        else:
                            reverse_rate_matching_parameters = list(reverse_rate_matching_parameters)

                        reverse_rate_constant = reverse_rate_matching_parameters[0]
                        reverse_rate_constant_value = parameters_values[reverse_rate_constant]
                        
                        if len( reverse_rate_matching_parameters ) == 2:

                            k_minus_1 = reverse_rate_matching_parameters[0]

                            k_minus_2 = reverse_rate_matching_parameters[1]

                            reverse_pattern = rf"(?:{re.escape(k_minus_1)}\s*([\+\-\*/])\s*{re.escape(k_minus_2)}|{re.escape(k_minus_2)}\s*([\+\-\*/])\s*{re.escape(k_minus_1)})"

                            reverse_match = re.search( reverse_pattern,  str(reverse_rate_expr) )

                            reverse_match_str = reverse_match.group()

                            if "*" in reverse_match_str:
                                reverse_rate_constant = k_minus_1 + " * " + k_minus_2
                                reverse_rate_constant_value = parameters_values[k_minus_1] *  parameters_values[k_minus_2]

                            elif "+" in reverse_match_str:
                                reverse_rate_constant = k_minus_1 + " + " + k_minus_2
                                reverse_rate_constant_value = parameters_values[k_minus_1] +  parameters_values[k_minus_2]

                            elif "-" in reverse_match_str:
                                reverse_rate_constant = k_minus_1 + " - " + k_minus_2
                                reverse_rate_constant_value = parameters_values[k_minus_1] -  parameters_values[k_minus_2]

                            elif "/" in reverse_match_str:
                                reverse_rate_constant = k_minus_1 + " / " + k_minus_2
                                reverse_rate_constant_value = parameters_values[k_minus_1] /  parameters_values[k_minus_2]

                            else:
                                raise ValueError(f"Operand not found in the reaction rate composition expression: {reverse_rate_expr} in line {inspect.currentframe().f_lineno()} in \"matrix_constructor.py\"")
                            
                        elif len( reverse_rate_matching_parameters ) > 2:

                            raise ValueError("Number of reaction rate constants are more than two which is not supported!")
                            


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

        kinetic_vector_length = 2 * reactions_number

        vector_of_kinetic_constants = np.zeros(kinetic_vector_length)

        for biomodel_reaction in biomodel_reactions:

            if biomodel_reaction.index != None:

                index = biomodel_reaction.index

                k_plus_value = biomodel_reaction.kinetic_forward_rate_constant_value if biomodel_reaction.kinetic_forward_rate_constant_value is not None else 0.0

                k_minus_value = biomodel_reaction.kinetic_reverse_rate_constant_value if biomodel_reaction.kinetic_reverse_rate_constant_value is not None else 0.0

                vector_of_kinetic_constants[index] = k_plus_value

                vector_of_kinetic_constants[index + reactions_number] = k_minus_value

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

        kinetic_rates_vector = self.kinetic_constants_vector_constructor(biomodel)

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

        reactions_number = Reaction.getCurrentIndex()

        ones_array = np.concatenate((np.ones(reactions_number), -np.ones(reactions_number)))

        result = (ones_array @ logn_kinetic_rates_vector).item()

        if abs(result) < 1e-6:

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are compatible with thermodynamic constraints", text_color="green")

            return True
        
        else:

            if printing.lower() == "on":
                utility.printer("\nCompatibility Check: ","The kinetic reaction rate constants are NOT compatible with thermodynamic constraints", text_color="red")
            
            return False