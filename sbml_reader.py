import libsbml
import utility
import re
from constants import *
from collections import defaultdict
import itertools
import sympy as sp
from classes.cBioMLReaction import *
from classes.cBioMLModel import *
from classes.cBioMLSpecies import *
from classes.cBioMLParameter import *
from classes.cBioMLSpeciesReference import *
from classes.cBioMLFunctionDefinition import *

import os
import exceptions
import time
import warnings

import model_checker




class SbmlReader:





    # ********************************
    # *           Function           *
    # ********************************
    def read_file(self, file_path):

        """
        Reads an SBML file using libSBML.

        :return: SBML model if successful, None otherwise.
        """

        self._file_name = os.path.basename(file_path)

        reader = libsbml.SBMLReader()
        document = reader.readSBML(file_path)
        if document.getNumErrors() > 0:
            utility.error_printer(f"The SBML file \"{self._file_name}\" contains {document.getNumErrors()} error(s).")
            utility.message_printer("\n>>>>> Model not read <<<<<<", color="red", style='bold')
            return None
        else:
            sbmodel = document.getModel()
            biomlmodel = BioMLModel(sbmodel.getId())
            biomlmodel.function_definitions = self._sbml_to_biomlmodel_function_definition_transfer(sbmodel)
            biomlmodel.species = self._sbml_to_biomlmodel_species_tranfer(sbmodel)
            biomlmodel.reactions = self._sbml_to_biomlmodel_reaction_tranfer(sbmodel, biomlmodel.function_definitions)
            biomlmodel.parameters = self._sbml_to_biomlmodel_parameter_transfer(sbmodel)
            biomlmodel.compartments = self._sbml_to_biomlmodel_compartments_transfer(sbmodel)

            _model_checker = model_checker.ModelChecker()

            if (_model_checker.check_mass_action_kinetics(biomlmodel, immediate_return=True)):

                try:

                    biomlmodel.is_mass_action = True

                    self._forward_reverse_rate_finder(biomlmodel)
                    
                    return biomlmodel
                
                except Exception:

                    return biomlmodel

            
            else:

                biomlmodel.is_mass_action = False
                
                return biomlmodel
        




    # ********************************
    # *           Function           *
    # ********************************
    def _sbml_to_biomlmodel_species_tranfer(self, libsbml_model):
        '''
        This function gets a SBML model, reads the required information for the species and creates a Species class for each one
        Then, it returns a list that contains the classes of species for this tool
        '''

        self._biomlmodel_species_list = []

        list_of_libsbml_species = libsbml_model.getListOfSpecies()

        for libsbml_species_class in list_of_libsbml_species:

            species_id = libsbml_species_class.getId()

            biomlmodel_species = BioMLSpecies(species_id)

            biomlmodel_species.initial_concentration = libsbml_species_class.getInitialConcentration()

            biomlmodel_species.compartment = libsbml_species_class.getCompartment()

            biomlmodel_species.charge = libsbml_species_class.getCharge()

            self._biomlmodel_species_list.append(biomlmodel_species)

        if not self._biomlmodel_species_list:
            utility.warning_printer("No species imported from SBML model!")

        return self._biomlmodel_species_list
    





    # ********************************
    # *           Function           *
    # ********************************
    def _sbml_to_biomlmodel_parameter_transfer(self, libsbml_model):

        biomlmodel_parameters_list = []

        libsbml_parameters_list = libsbml_model.getListOfParameters()

        for libsbml_parameter_class in libsbml_parameters_list:

            parameter_id = libsbml_parameter_class.getId()

            parameter_value = libsbml_parameter_class.getValue()

            biomlmodel_parameter = BioMLParameter(parameter_id)

            biomlmodel_parameter.value = parameter_value

            biomlmodel_parameters_list.append(biomlmodel_parameter)

        try:

            if not biomlmodel_parameters_list:
                utility.warning_printer("No GLOBAL parameters found in the SBML model!")
                time.sleep(1)

                biomlmodel_parameters_list = self._sbml_local_parameter_finder(libsbml_model)  
            
        except exceptions.LocalParameterConflict as e:

            utility.warning_printer(e)

            return biomlmodel_parameters_list
        
        except exceptions.EmptyList as e:

            utility.warning_printer(e)

            return biomlmodel_parameters_list


        return biomlmodel_parameters_list
    





    # ********************************
    # *           Function           *
    # ********************************
    def _sbml_to_biomlmodel_reaction_tranfer(self, libsbml_model, function_definitions):

        biomlmodel_reactions_list = []

        libsbml_reactions = libsbml_model.getListOfReactions()

        for libsbml_reaction_class in libsbml_reactions:

            biomlmodel_products_list =[]
            biomlmodel_reactants_list = []

            libsbml_reactants = libsbml_reaction_class.getListOfReactants()

            libsbml_products = libsbml_reaction_class.getListOfProducts()

            reaction_id = libsbml_reaction_class.getId()

            index = BioMLReaction.get_current_index()

            biomlmodel_reaction = BioMLReaction(reaction_id)

            biomlmodel_reaction.reversible = libsbml_reaction_class.getReversible()

            libsbml_klaw = libsbml_reaction_class.getKineticLaw()

            if libsbml_klaw:

                try:

                    biomlmodel_reaction.kinetic_law = libsbml_klaw.getFormula()

                except Exception as e:
                    raise ValueError(f"The formula can't be read from the SBML reaction law. The message from libsbl is: " + str(e))

            else:

                raise ValueError(f"Reaction {reaction_id} does not have a reaction rate formula")

            biomlmodel_reaction.expanded_kinetic_law , _ = SbmlReader._expand_formula(biomlmodel_reaction.kinetic_law, function_definitions)

            sbml_level = libsbml_model.getLevel()

            if sbml_level == 3:

                sbml_local_parameters = libsbml_reaction_class.getKineticLaw().getListOfLocalParameters()

                local_parameters = []

                for sbml_local_parameter in sbml_local_parameters:

                    local_parameter_id = sbml_local_parameter.getId()

                    local_parameter_value = sbml_local_parameter.getValue()

                    biomlmodel_parameter = BioMLParameter(local_parameter_id)

                    biomlmodel_parameter.value = local_parameter_value

                    local_parameters.append(biomlmodel_parameter)

                biomlmodel_reaction.local_parameters = local_parameters

            elif sbml_level == 1 or sbml_level == 2:

                sbml_local_parameters = libsbml_reaction_class.getKineticLaw().getListOfParameters()

                local_parameters = []

                for sbml_local_parameter in sbml_local_parameters:

                    local_parameter_id = sbml_local_parameter.getId()

                    local_parameter_value = sbml_local_parameter.getValue()

                    biomlmodel_parameter = BioMLParameter(local_parameter_id)

                    biomlmodel_parameter.value = local_parameter_value

                    local_parameters.append(biomlmodel_parameter)

                biomlmodel_reaction.local_parameters = local_parameters

            for libsbml_reactant_class in libsbml_reactants:

                id = libsbml_reactant_class.getSpecies()

                if id == "empty":

                    BioMLReaction.reset_counter(index)

                    biomlmodel_reaction.reset_index()

                    biomlmodel_reaction.boundary_condition = True

                for biomlmodel_species in self._biomlmodel_species_list:

                    if id == biomlmodel_species.ID:

                        biomlmodel_species_reference = BioMLSpeciesReference(biomlmodel_species)

                        biomlmodel_species_reference.reaction_id = reaction_id

                        biomlmodel_species_reference.stoichiometry = libsbml_reactant_class.getStoichiometry()

                        biomlmodel_reactants_list.append(biomlmodel_species_reference)

            for libsbml_product_class in libsbml_products:

                id = libsbml_product_class.getSpecies()

                if id == "empty":

                    BioMLReaction.reset_counter(index)

                    biomlmodel_reaction.reset_index()

                    biomlmodel_reaction.boundary_condition = True

                for biomlmodel_species in self._biomlmodel_species_list:

                    if id == biomlmodel_species.ID:

                        biomlmodel_species_reference = BioMLSpeciesReference(biomlmodel_species)

                        biomlmodel_species_reference.reaction_id = reaction_id

                        biomlmodel_species_reference.stoichiometry = libsbml_product_class.getStoichiometry()

                        biomlmodel_products_list.append(biomlmodel_species_reference)

            biomlmodel_reaction.reactants = biomlmodel_reactants_list

            biomlmodel_reaction.products = biomlmodel_products_list

            biomlmodel_reactions_list.append(biomlmodel_reaction)

        if not biomlmodel_reactions_list:
           utility.warning_printer("There are no reactions defined in the SBML file!!!")

        return biomlmodel_reactions_list
    





    # ********************************
    # *           Function           *
    # ********************************
    def _sbml_to_biomlmodel_function_definition_transfer(self, libsbml_model):

        sbml_function_definitons = [libsbml_model.getFunctionDefinition(i)
                                    for i in range(libsbml_model.getNumFunctionDefinitions())]
        
        function_definitions = [BioMLFunctionDefinition(sb)
                                for sb in sbml_function_definitons]
        
        return function_definitions
    






    # ********************************
    # *           Function           *
    # ********************************
    def _sbml_local_parameter_finder(self, libsbml_model):

        local_parameters_strings = []

        libsbml_reactions_list = libsbml_model.getListOfReactions()

        conflicting_parameter_IDs = False

        for libsbml_reaction_class in libsbml_reactions_list:

            libsbml_reaction_parameters_list = libsbml_reaction_class.getKineticLaw().getListOfParameters()

            for libsbml_reaction_parameter_class in libsbml_reaction_parameters_list:

                parameter_id = libsbml_reaction_parameter_class.getId()

                if parameter_id not in local_parameters_strings:

                    local_parameters_strings.append(parameter_id)

                else:

                    conflicting_parameter_IDs = True

        if conflicting_parameter_IDs is False:
                    
            local_biomlmodel_parameters_list = []

            for libsbml_reaction_class in libsbml_reactions_list:

                libsbml_reaction_parameters_list = libsbml_reaction_class.getKineticLaw().getListOfParameters()
                    
                for libsbml_reaction_parameter_class in libsbml_reaction_parameters_list:
                        
                    parameter_id = libsbml_reaction_parameter_class.getId()

                    parameter_value = libsbml_reaction_parameter_class.getValue()

                    biomlmodel_parameter = BioMLParameter(parameter_id)

                    biomlmodel_parameter.value = parameter_value

                    local_biomlmodel_parameters_list.append(biomlmodel_parameter)

            if not local_biomlmodel_parameters_list:
                raise exceptions.EmptyList("There are no local parameters!")
            else:
                utility.message_printer("\nLocal parameters are stored as global parameters, too!", color='blue')
                time.sleep(1)

            return local_biomlmodel_parameters_list

        else:

            raise exceptions.LocalParameterConflict("Global parameters list cannot be created from local parameters!\n         There are conflictions in local parameter names in different reactions!")
        






    # ********************************
    # *           Function           *
    # ********************************
    def _sbml_to_biomlmodel_compartments_transfer(self, libsbml_model):

        compartments = [comp.getId()
                        for comp in libsbml_model.getListOfCompartments()]
        
        return compartments
    




    # ********************************
    # *           Function           *
    # ********************************
    def _forward_reverse_rate_finder(self, biomlmodel, printing = "off") -> dict:
        """
        Input of this function is a biomlmodel class instance
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
        def _get_forward_reverse_rate_expressions(expression):

            # Separate positive and negative terms manually
            rates = defaultdict(list)

            # rates = {
            #     "forward_rate": [],
            #     "reverse_rate": []
            # }

            # Check for invalid or empty expressions
            if expression is None or expression == 0:
                return rates  # Nothing to do

            if not isinstance(expression, sp.Expr):
                raise TypeError("Expected a sympy expression.")

            if expression.is_Mul:
                flag = True
                for argument in expression.args:
                    if argument.is_Add:
                        flag = False
                        for arg in expression.args:
                            if "-" in str(argument):
                                rates["reverse_rate"].append(arg)
                            else:
                                rates["forward_rate"].append(arg)
                if flag:
                    rates["forward_rate"].append(expression)

            elif expression.is_Add:
                for argument in expression.args:
                    if "-" in str(argument):
                        rates["reverse_rate"].append(argument)
                    else:
                        rates["forward_rate"].append(argument)

            else:
                rates["forward_rate"].append(expression)

            return rates
        # --------------------------------------------------
        # --------------------------------------------------

        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def _power_operator_finder(expression, parameters, parameters_values, new_parameters=None, reaction_rate_constant='', reaction_rate_constant_value=None, number_of_recursions = 0):

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

                    used_elements.add(pair[0])
                    
                    used_elements.add(pair[1])

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

            return _power_operator_finder(expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value, number_of_recursions+1)
        # --------------------------------------------------
        # --------------------------------------------------


        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def _parameters_finder(expression, parameters, parameters_values, new_parameters = None, reaction_rate_constant = "", reaction_rate_constant_value=None, number_of_recursions=0):

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

                    used_elements.add(pair[0])

                    used_elements.add(pair[1])

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

            return _parameters_finder(expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value, number_of_recursions+1)
        # --------------------------------------------------
        # --------------------------------------------------


        empty_global_parameters = False     # a flag to indicate if there are no global parameters

        if biomlmodel is None:
            raise exceptions.NoModel("No biomlmodel has been read!!!")


        species_classes_list = biomlmodel.get_list_of_species()
        parameter_classes_list = biomlmodel.get_list_of_parameters()
        reaction_classes_list = biomlmodel.get_list_of_reactions()
        function_definitions_list = biomlmodel.get_list_of_function_definitions()
        compartments_list = biomlmodel.get_list_of_compartments()

        empty_species_list = False

        if not species_classes_list:
            empty_species_list = True

        if not parameter_classes_list:
            empty_global_parameters = True

        if not reaction_classes_list:
            raise ValueError("Reaction list is empty.")



        string_to_sympy_symbols = {} # THERE ARE SOME BUILT-IN FUNCTIONS IN SYMPY THAT INTERFERE WITH VARIABLE NAMES AND THE STRING TO ... \
                                        # SYMPY EXPRESSIONS AND THE EXECUTION FAILS. TO PREVENT IT, I HAVE DEFINED SOME BUILT-IN FUNCTIONS AS ... \
                                        # SYMPY VARIABLES BEFOREHAND

        for individual_species_class in species_classes_list:
            species_name = individual_species_class.get_id()
            string_to_sympy_symbols[species_name] = sp.symbols(species_name)

            string_to_sympy_symbols.update({
                name: sp.symbols(name) for name in compartments_list
            })


        parameters_values = {}  # This list stores the names of the parameters and their values

        for individual_parameter_class in parameter_classes_list:

            parameter_name = individual_parameter_class.get_id()
            parameter_value = individual_parameter_class.get_value()
            string_to_sympy_symbols[parameter_name] = sp.symbols(parameter_name)
            parameters_values[parameter_name] = parameter_value


        for individual_reaction_class in reaction_classes_list:

            if individual_reaction_class.local_parameters is not None:

                local_parameters_values = {}

                local_parameters_list = individual_reaction_class.local_parameters

                for local_parameter_class in local_parameters_list:

                    local_parameter_name = local_parameter_class.get_id()
                    local_parameter_value = local_parameter_class.get_value()
                    string_to_sympy_symbols[local_parameter_name] = sp.symbols(local_parameter_name)
                    local_parameters_values[local_parameter_name] = local_parameter_value

                parameters_values.update(local_parameters_values)
            
            else:

                if empty_global_parameters == True:
                    raise ValueError(f"There are no global parameters defined for the model, nor are there any local parameters defined for {reaction_name}!")
                
            reactant_classes_list = individual_reaction_class.get_list_of_reactants()
                
            for individual_reactant_class in reactant_classes_list:
                reactant_name = individual_reactant_class.get_id()
                string_to_sympy_symbols[reactant_name] = sp.symbols(reactant_name)
                empty_species_list = False

            product_classes_list = individual_reaction_class.get_list_of_products()
                
            for individual_product_class in product_classes_list:
                product_name = individual_product_class.get_id()
                string_to_sympy_symbols[product_name] = sp.symbols(product_name)
                empty_species_list = False

            if empty_species_list:
                raise ValueError(f"There are no species to be checked for reaction {reactant_name}")


            reaction_name = individual_reaction_class.get_id()
            reaction_rate_formula = individual_reaction_class.get_kinetic_law()

            reaction_rate_formula, function_symbols = SbmlReader._expand_formula(reaction_rate_formula, function_definitions_list)

            string_to_sympy_symbols.update(function_symbols)

            sp_reaction_rate_formula = sp.sympify(reaction_rate_formula, locals = string_to_sympy_symbols, evaluate = False)

            simplified_formula = sp.cancel(sp_reaction_rate_formula)

            expanded_formula = sp.expand(simplified_formula)

            if printing.lower() == 'on':
                utility.printer(f"\nThe simplified reaction rate expression for reaction {reaction_name} is:\n", simplified_formula)

            

            forward_reverse_rate_equations = _get_forward_reverse_rate_expressions(expanded_formula)

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

            if not forward_rate_matching_parameters:
                raise ValueError(f"There is not a forward reaction rate constant in reaction {reaction_name} with expanded kinetic law {str(expanded_formula)}")

            forward_rate_constant = forward_rate_matching_parameters[0]
            forward_rate_constant_value = parameters_values[forward_rate_constant]

            if len( forward_rate_matching_parameters ) > 2 and len( forward_rate_matching_parameters ) < MAX_RECURSION:

                parameters_list = [str(p) for p in forward_rate_matching_parameters]

                expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value = _power_operator_finder(str(forward_rate_expression), parameters_list, parameters_values)

                expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value = _parameters_finder(expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value)

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

                if not forward_rate_matching_parameters:
                    raise ValueError(f"There is not a forward reaction rate constant in reaction {reaction_name} with expanded kinetic law: {str(expanded_formula)}")

                temp_forward_rate_constant = forward_rate_matching_parameters[0]
                temp_forward_rate_constant_value = parameters_values[temp_forward_rate_constant]

                if len( forward_rate_matching_parameters ) > 2 and len( forward_rate_matching_parameters ) <= MAX_RECURSION:

                    parameters_list = [str(p) for p in forward_rate_matching_parameters]

                    expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value = _power_operator_finder(str(forward_rate_expression), parameters_list, parameters_values)

                    expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value = _parameters_finder(expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value)

                    message = f"\nThe forward kinetic rate constant for reaction {reaction_name} has more than one variable: {forward_rate_constant}"

                    utility.add_warning(message)

                elif len( forward_rate_matching_parameters ) > MAX_RECURSION:

                    raise ValueError(f"Number of reaction rate constants are more than {MAX_RECURSION} which is not supported!")
                
                forward_rate_constant += " + " + temp_forward_rate_constant
                forward_rate_constant_value += temp_forward_rate_constant_value

                message = f"\nThe forward kinetic rate constant for reaction {reaction_name} has more than one variable: {forward_rate_constant}"

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

                    if not forward_rate_matching_parameters:
                        raise ValueError(f"There is not a reverse reaction rate constant in reaction {reaction_name}")

                    reverse_rate_constant = reverse_rate_matching_parameters[0]
                    reverse_rate_constant_value = parameters_values[reverse_rate_constant]

                    if len( reverse_rate_matching_parameters ) > 2 and len( reverse_rate_matching_parameters ) <= MAX_RECURSION:

                        parameters_list = [str(p) for p in reverse_rate_matching_parameters]

                        expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value = _power_operator_finder(str(reverse_rate_expression), parameters_list, parameters_values)

                        expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value = _parameters_finder(expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value)

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

                        if not forward_rate_matching_parameters:
                            raise ValueError(f"There is not a reverse reaction rate constant in reaction{reaction_name}")

                        temp_reverse_rate_constant = reverse_rate_matching_parameters[0]
                        temp_reverse_rate_constant_value = parameters_values[temp_reverse_rate_constant]
                        
                        
                        if len( reverse_rate_matching_parameters ) > 2 and len( reverse_rate_matching_parameters ) <= MAX_RECURSION:

                            parameters_list = [str(p) for p in reverse_rate_matching_parameters]

                            expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value = _power_operator_finder(str(reverse_rate_expression), parameters_list, parameters_values)

                            expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value = _parameters_finder(expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value)

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
    







    @staticmethod
    def _expand_formula(formula, function_definitions,
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
                    if arg:
                        input_symbols_dict[arg] = sp.symbols(arg)
                function_formula = function_definition.formula
                for formal_arg, call_arg in zip(function_definition.arguments, input_arguments):
                    function_formula = re.sub(rf'\b{re.escape(formal_arg)}\b', call_arg, function_formula)
                formula = formula.replace(function_call, function_formula)
        if not done:
            return SbmlReader._expand_formula(formula, function_definitions,
                        num_recursions=num_recursions+1)
        return formula, input_symbols_dict