import numpy as np
import exceptions
import os
import sympy as sp
from classes.cReaction import *
from classes.cSpecies import *
import utility
import random
import warnings

    



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

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            utility.error_printer("\nError: ", e)
            print("Unable to complete the query\!")
            return

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

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            utility.error_printer("\nError: ", e)
            print("Unable to complete the query\!")
            return

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

        try:

            if len(species_list) == 0:
                raise exceptions.EmptyList("There are no species in this model.")
            
            if len(reactions_list) == 0:
                raise exceptions.EmptyList("There are no reactions in this model.")
            
        except exceptions.EmptyList as e:
            utility.error_printer("\nError: ", e)
            print("Unable to complete the query\!")
            return

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

        if not species_classes_list:
            raise ValueError("Species list is empty.")

        if not parameter_classes_list:
            raise ValueError("Parameter list is empty.")

        if not reaction_classes_list:
            raise ValueError("Reaction list is empty.")



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

                utility.error_printer("\nSympify Error: ", e)
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")

            except ValueError as v:
                utility.error_printer("\nValue Error: ", e)
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")

            except Exception as e:
                utility.error_printer("\nUnexpected Error: ", e)
                print(f"\nEquation couldn't be converted to Sympy expression for reaction: {reaction_name}")
                print(f"\nThe string for equation is: {reaction_rate_formula}")

            else:

                if printing.lower() == 'on':
                    utility.printer("\nThe reaction rate expression is:\n", sp_reaction_rate_formula, text_color="yellow")

            simplified_formula = sp.simplify(sp_reaction_rate_formula)

            forward_reverse_rate_equations = get_forward_reverse_rate_expressions(simplified_formula)

            forward_variables_symbols = sp.sympify(forward_reverse_rate_equations.get("forward_rate"), locals = string_to_sympy_symbols, evaluate = False).free_symbols
            forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
            forward_rate_constant = next(iter(set(forward_variables_as_strings) & set(parameters_values.keys())), None)

            if forward_rate_constant:

                try:

                    individual_reaction_class.kinetic_forward_rate_constant = forward_rate_constant
                    individual_reaction_class.kinetic_forward_rate_constant_value = parameters_values[str(forward_rate_constant)]

                except ValueError as e:
                    utility.error_printer(f"\nERROR: ", e)
                    utility.message_printer(f"Forward reaction rate constant can't be found for {reaction_name}", style="normal")

                reverse_rate_constant = None

                if forward_reverse_rate_equations.get("reverse_rate"):

                    reverse_variables_symbols = sp.sympify(forward_reverse_rate_equations.get("reverse_rate"), locals = string_to_sympy_symbols, evaluate = False).free_symbols
                
                    reverse_variables_as_strings = [str(symbol) for symbol in reverse_variables_symbols]

                    reverse_rate_constant = next(iter(set(reverse_variables_as_strings) & set(parameters_values.keys())), None)

                    try:

                        individual_reaction_class.kinetic_reverse_rate_constant = reverse_rate_constant
                        individual_reaction_class.kinetic_reverse_rate_constant_value = parameters_values[str(reverse_rate_constant)]

                    except ValueError as e:
                        utility.error_printer(f"\nERROR: ", e)
                        utility.message_printer(f"Reverse reaction rate constant can't be found for {reaction_name}", style="normal")

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

            else:

                if forward_reverse_rate_equations.get("reverse_rate"):

                    forward_variables_symbols = sp.sympify(simplified_formula, locals = string_to_sympy_symbols, evaluate = False).free_symbols
                    forward_variables_as_strings = [str(symbol) for symbol in forward_variables_symbols]
                    common_rate_constant = next(iter(set(forward_variables_as_strings) & set(parameters_values.keys())), None)

                    forward_rate_constant = common_rate_constant

                    reverse_rate_constant = common_rate_constant

                    try:

                        individual_reaction_class.kinetic_forward_rate_constant = forward_rate_constant
                        individual_reaction_class.kinetic_forward_rate_constant_value = parameters_values[str(forward_rate_constant)]

                        individual_reaction_class.kinetic_reverse_rate_constant = reverse_rate_constant
                        individual_reaction_class.kinetic_reverse_rate_constant_value = parameters_values[str(reverse_rate_constant)]

                    except ValueError as e:
                        utility.error_printer(f"\nERROR: ", e)
                        utility.message_printer(f"Reaction rate constant can't be found for {reaction_name}", style="normal")

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
        
            self.forward_reverse_rate_finder(biomodel)

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

                return False

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