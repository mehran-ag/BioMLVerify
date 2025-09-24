import libsbml
import _modules._utility as utility
import re
from _modules._constants import *
from collections import defaultdict
import itertools
import sympy as sp
from _classes.cBioMLReaction import *
from _classes.cBioMLModel import *
from _classes.cBioMLSpecies import *
from _classes.cBioMLParameter import *
from _classes.cBioMLSpeciesReference import *
from _classes.cBioMLFunctionDefinition import *

import os
import _modules._exceptions as exceptions
import time

import chemparse as chp
import libchebipy as chb

import math

import _modules._model_checker as model_checker




class SbmlReader:





    # ********************************
    # *           Function           *
    # ********************************
    def read_file(self, file_path: str) -> BioMLModel:

        """
            Reads an SBML file using libSBML.

            Args:
                file_path (str): the full path to the file including file name

            Returns:
                BioMLModel: An instance of BioML Model where all contents of SBML file have been converted to BioML specific counterparts
        """

        self._file_name = os.path.basename(file_path)

        reader = libsbml.SBMLReader()
        document = reader.readSBML(file_path)
        if document.getNumErrors() > 0:
            utility.error_printer(f"The SBML file \"{self._file_name}\" contains {document.getNumErrors()} error(s).")

            for i in range(document.getNumErrors()):
                error = document.getError(i)
                utility.error_printer(f"{error.getMessage()}")


            utility.message_printer("\n>>>>> Model not read <<<<<<", color="red", style='bold')
            return None
        else:
            sbmodel = document.getModel()
            biomlmodel = BioMLModel(sbmodel.getId())
            biomlmodel.function_definitions = self._transfer_sbml_function_definitions_to_biomlmodel(sbmodel)
            biomlmodel.species = self._transfer_sbml_species_to_biomlmodel(sbmodel)
            biomlmodel.reactions = self._transfer_sbml_reactions_to_biomlmodel(sbmodel, biomlmodel.function_definitions)
            biomlmodel.parameters = self._transfer_sbml_parameters_to_biomlmodel(sbmodel)
            biomlmodel.compartments = self._get_list_of_sbml_compartments(sbmodel)

            _model_checker = model_checker.ModelChecker()

            if (_model_checker.check_mass_action_kinetics(biomlmodel, immediate_return=True)):

                try:

                    biomlmodel.is_mass_action = True

                    self._find_forward_reverse_rate_constants(biomlmodel)
                    
                    return biomlmodel
                
                except Exception:

                    return biomlmodel

            
            else:

                biomlmodel.is_mass_action = False
                
                return biomlmodel
        




    # ********************************
    # *           Function           *
    # ********************************
    def _transfer_sbml_species_to_biomlmodel(self, libsbml_model: libsbml.Model) -> list[BioMLSpecies]:
        """
            Converts species from an SBML model into BioMLSpecies instances.

            Parses the input `libsbml` model to extract species information and creates a `BioMLSpecies` instance for each species.
            Returns a list of all such instances for use in this tool.

            Args:
                libsbml_model (libsbml.Model): An SBML model from the `libsbml` package.

            Returns:
                list[object]: A list of BioMLSpecies instances created from the SBML species.
        """


        self._biomlmodel_species_list = []

        list_of_libsbml_species = libsbml_model.getListOfSpecies()

        for libsbml_species_class in list_of_libsbml_species:

            species_id = libsbml_species_class.getId()

            biomlmodel_species = BioMLSpecies(species_id)

            biomlmodel_species.initial_concentration = libsbml_species_class.getInitialConcentration()

            biomlmodel_species.compartment = libsbml_species_class.getCompartment()

            annotations = SbmlReader._get_chebi_annotations(libsbml_species_class)

            if annotations:

                biomlmodel_species.annotations['chebi'] = annotations

                for annotation in annotations:

                    formula, charge, composition = SbmlReader._parse_using_chebi(annotation)

                    if formula is not None or charge is not None or composition is not None:
                        break

                if formula is not None:
                    biomlmodel_species.compound = formula

                if composition is not None:
                    biomlmodel_species.composition = composition

                if charge is not None:

                    if not math.isnan(charge):

                        if libsbml_species_class.getCharge():
                            biomlmodel_species.charge = libsbml_species_class.getCharge()
                        else:
                            if charge:
                                biomlmodel_species.charge = charge

                    else:

                        biomlmodel_species.charge = 0

            else:

                sbml_metaid = libsbml_species_class.getMetaId()

                if sbml_metaid.split('_')[0] == 'va':

                    compound_part = sbml_metaid.split('_')[1]

                    if compound_part[0] == '-':
                        charge_sign = -1
                        compound_part = compound_part[1:]

                    else:
                        charge_sign = +1

                    compound_code = compound_part.split('-')[0]

                    sbml_species_compound, sbml_species_charge = SbmlReader._parse_compound_code(compound_code)

                    biomlmodel_species.charge = charge_sign * sbml_species_charge

                    if len( compound_part.split('-') ) > 1:

                        sbml_species_comp_code = compound_part.split('-')[1]

                        sbml_species_composition = self._parse_molecule_units(sbml_species_comp_code)

                    else:

                        sbml_species_composition = {sbml_species_compound: 1}

                    biomlmodel_species.compound = sbml_species_compound

                    biomlmodel_species.composition = sbml_species_composition


            self._biomlmodel_species_list.append(biomlmodel_species)

        if not self._biomlmodel_species_list:
            utility.warning_printer("No species imported from SBML model!")

        return self._biomlmodel_species_list




    # ********************************
    # *           Function           *
    # ********************************
    def _transfer_sbml_parameters_to_biomlmodel(self, libsbml_model: libsbml.Model) -> list[BioMLParameter]:
        """
            Converts parameters from an SBML model into BioMLParameter instances.

            Parses the input `libsbml` model to extract parameter information and creates a `BioMLParameter` instance for each parameter.
            Returns a list of all such instances for use in this tool.

            Args:
                libsbml_model (libsbml.Model): An SBML model from the `libsbml` package.

            Returns:
                list[object]: A list of BioMLParameter instances created from the SBML parameters.
        """

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

                biomlmodel_parameters_list = self._find_sbml_local_parameters(libsbml_model)  
            
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
    def _transfer_sbml_reactions_to_biomlmodel(self, libsbml_model: libsbml.Model, bioml_function_definitions: list[BioMLFunctionDefinition]) -> list[BioMLReaction]:
        """
            Converts reactions from an SBML model into BioMLReaction instances.

            Parses the input `libsbml` model to extract reaction information and creates a `BioMLReaction` instance for each reaction.
            Returns a list of all such instances for use in this tool.

            Args:
                libsbml_model (libsbml.Model): An SBML model from the `libsbml` package.

            Returns:
                list[object]: A list of BioMLReaction instances created from the SBML reactions.
        """

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

            biomlmodel_reaction.expanded_kinetic_law , _ = SbmlReader._expand_formula(biomlmodel_reaction.kinetic_law, bioml_function_definitions)

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
    def _transfer_sbml_function_definitions_to_biomlmodel(self, libsbml_model: libsbml.Model) -> list[BioMLFunctionDefinition]:
        """
            Converts function definitions from an SBML model into BioMLFunctionDefinition instances.

            Parses the input `libsbml` model to extract function definition information and creates a `BioMLFunctionDefinition` instance for each function definition.
            Returns a list of all such instances for use in this tool.

            Args:
                libsbml_model (libsbml.Model): An SBML model from the `libsbml` package.

            Returns:
                list[object]: A list of BioMLFunctionDefinition instances created from the SBML function definitions.
        """

        sbml_function_definitons = [libsbml_model.getFunctionDefinition(i)
                                    for i in range(libsbml_model.getNumFunctionDefinitions())]
        
        bioml_function_definitions = [BioMLFunctionDefinition(sb)
                                for sb in sbml_function_definitons]
        
        return bioml_function_definitions
    






    # ********************************
    # *           Function           *
    # ********************************
    def _find_sbml_local_parameters(self, libsbml_model: libsbml.Model) -> list[BioMLParameter]:
        """
            Converts local parameters from an SBML model into BioMLParameter instances.

            Parses the input `libsbml` model to extract parameter information and creates a `BioMLParameter` instance for each local parameter.
            Returns a list of all such instances for use in this tool.

            Args:
                libsbml_model (libsbml.Model): An SBML model from the `libsbml` package.

            Returns:
                list: A list of BioMLParameter instances created from the SBML local parameters.
        """

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
    def _get_list_of_sbml_compartments(self, libsbml_model: libsbml.Model) -> list[str]:

        compartments = [comp.getId()
                        for comp in libsbml_model.getListOfCompartments()]
        
        return compartments
    




    # ********************************
    # *           Function           *
    # ********************************
    def _find_forward_reverse_rate_constants(self, biomlmodel: BioMLModel, printing: bool = False) -> dict:
        """
            Extracts the rate constants for each reaction in the given model, 
            including both forward and reverse directions.
            Then, "kinetic_forward_rate_constant" and "kinetic_reverse_rate_constant" variables of each reaction and their values are updated.

            Each reaction's rate constants are stored in a dictionary with 'forward' and 'reverse' 
            as keys and their corresponding values as floats.

            These per-reaction dictionaries are then stored in a master dictionary, 
            `reaction_to_rate_constants`, where the keys are reaction names and the values 
            are the corresponding rate constant dictionaries.

            The final dictionary is returned.

            Args:
                biomlmodel (BioMLModel): A model of the BioML class containing species and reactions.
                printing (bool): If True, displays the reaction rate constants for each reaction as they are found

            Returns:
                dict: A dictionary containing reaction names as keys and another dictionary containing forard and reverse reaction rate names mapped to their values as the dictionary's value
        """

        reaction_to_rate_constants = {} # This dictionary maps reaction names (IDs) to a dictionary of forward and reverse rate constants and their values



        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def _get_forward_reverse_rate_expressions(expression):
            """
                Classifies terms in a SymPy expression into forward and reverse reaction rates.

                This function takes a symbolic expression typically representing a rate law
                or net reaction rate and separates its terms into forward and reverse contributions.
                Positive terms are interpreted as forward rates, while negative terms are
                considered reverse rates.

                Parameters:
                    expression (sp.Expr): A SymPy expression representing a rate or combination
                                        of rates. Can include addition and multiplication of terms.

                Returns:
                    dict: A dictionary with two keys:
                        - 'forward_rate': a list of SymPy expressions contributing positively.
                        - 'reverse_rate': a list of SymPy expressions contributing negatively.

                Notes:
                    - If the expression is zero or None, an empty dictionary is returned.
                    - If the expression is a multiplication of terms and one or more of them
                    is an addition (e.g., (A + B)*C), each additive term is classified individually.
                    - The check for negativity is performed by searching for a "-" in the string
                    representation of the term, which may not be robust for all SymPy expressions.
                    Consider improving this with SymPy's `.is_negative` if needed.
                
                Raises:
                    TypeError: If the input is not a SymPy expression.
            """

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
        def _find_power_operator(expression: str, parameters: list[str], parameters_values: dict, new_parameters: list[str] = None, reaction_rate_constant: str = '', reaction_rate_constant_value: float = None, number_of_recursions: int = 0):
            """
                Finds a parameters powered to a variable or constant and repalces it with a new paramater recursively until all parameters powered are found.

                Args:
                    expression (str): The equation to be looked up
                    parameters list[str]: a list of parameters in the reactions
                    parameters_values (dict): a dictionary mapping parameter names to their values
                    new_parameters (list[str]): the new parameters created to replace the powered parameter
                    reaction_rate_constant (str): this string stores the powered parameter and adds new ones to keep them as the paramater for the reation rate expression
                    reaction_rate_constant_value (float): the value of the created new parameter
                    number_of_recursions (int): this number is checked in every recursion to avoid infinite recursions

                Returns:
                    expression (str): The equation to be looked up
                    parameters list[str]: a list of parameters in the reactions
                    parameters_values (dict): a dictionary mapping parameter names to their values
                    new_parameters (list[str]): the new parameters created to replace the powered parameter
                    reaction_rate_constant (str): this string stores the powered parameter and adds new ones to keep them as the paramater for the reation rate expression
                    reaction_rate_constant_value (float): the value of the created new parameter
                    number_of_recursions (int): this number is checked in every recursion to avoid infinite recursions
            """

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

            return _find_power_operator(expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value, number_of_recursions+1)
        # --------------------------------------------------
        # --------------------------------------------------


        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # *      Internal Function       *
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        def _find_parameters(expression: str, parameters: list[str], parameters_values: dict, new_parameters: list = None, reaction_rate_constant: str = "", reaction_rate_constant_value: float = None, number_of_recursions: int = 0):
            """
                Finds two parameters added, multiplied, subtracted or divided and repalces them with a new paramater recursively until all parameters powered are found.

                Args:
                    expression (str): The equation to be looked up
                    parameters list[str]: a list of parameters in the reactions
                    parameters_values (dict): a dictionary mapping parameter names to their values
                    new_parameters (list[str]): the new parameters created to replace the powered parameter
                    reaction_rate_constant (str): this string stores the powered parameter and adds new ones to keep them as the paramater for the reation rate expression
                    reaction_rate_constant_value (float): the value of the created new parameter
                    number_of_recursions (int): this number is checked in every recursion to avoid infinite recursions

                Returns:
                    expression (str): The equation to be looked up
                    parameters list[str]: a list of parameters in the reactions
                    parameters_values (dict): a dictionary mapping parameter names to their values
                    new_parameters (list[str]): the new parameters created to replace the powered parameter
                    reaction_rate_constant (str): this string stores the powered parameter and adds new ones to keep them as the paramater for the reation rate expression
                    reaction_rate_constant_value (float): the value of the created new parameter
                    number_of_recursions (int): this number is checked in every recursion to avoid infinite recursions
            """

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

            return _find_parameters(expression, parameters, parameters_values, new_parameters, reaction_rate_constant, reaction_rate_constant_value, number_of_recursions+1)
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

            if printing:
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

                expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value = _find_power_operator(str(forward_rate_expression), parameters_list, parameters_values)

                expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value = _find_parameters(expression, parameters, parameters_values, new_parameters, forward_rate_constant, forward_rate_constant_value)

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

                    expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value = _find_power_operator(str(forward_rate_expression), parameters_list, parameters_values)

                    expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value = _find_parameters(expression, parameters, parameters_values, new_parameters, temp_forward_rate_constant, temp_forward_rate_constant_value)

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

                        expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value = _find_power_operator(str(reverse_rate_expression), parameters_list, parameters_values)

                        expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value = _find_parameters(expression, parameters, parameters_values, new_parameters, reverse_rate_constant, reverse_rate_constant_value)

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

                            expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value = _find_power_operator(str(reverse_rate_expression), parameters_list, parameters_values)

                            expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value = _find_parameters(expression, parameters, parameters_values, new_parameters, temp_reverse_rate_constant, temp_reverse_rate_constant_value)

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

                        if printing:

                            utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")

                        reaction_to_rate_constants[reaction_name] = {"forward": forward_rate_constant_value}

                    else:

                        if printing:

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

                            if printing:

                                utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")

                            reaction_to_rate_constants[reaction_name] = {"forward": parameters_values[str(forward_rate_constant)]}

                        else:

                            if printing:

                                utility.printer("\nForward rate constant is:", forward_rate_constant, text_style="bold")
                                utility.printer("Reverse rate constant is:", reverse_rate_constant, text_style="bold")


                            reaction_to_rate_constants[reaction_name] = {"forward": parameters_values[str(forward_rate_constant)], "reverse": parameters_values[str(reverse_rate_constant)]}



        return reaction_to_rate_constants
    







    @staticmethod
    def _expand_formula(formula: str, function_definitions: list[BioMLFunctionDefinition],
            num_recursions: int = 0) -> tuple[str, dict]:
        """
            Expands the kinetics formula, replacing a variable which is defiend by a function with their function definitions

            Args:
                formula (str): expansion of the kinetic law
                
                function_definitions (list[BioMLFunctionDefinition]): list of BioMLFunctionDefinition instances
                num_recursion (int): an integer counter to prevent infinite recursions
            
            Returns:
                str: the updated formula
                dict: a dictionary mapping the variables of the new functions to their Sympy Symbols
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
    


    @staticmethod
    def _get_chebi_annotations(libsbml_species: libsbml.Species):
        """
            Finds the ChEBI codes in the annotation node of the species and stores them in a list

            Args:
                species (libsbml.Species: A libsbml species instance

            Returns:
                list[str]: ChEBI codes (as digits only) stored in a list
        """

        annotation = libsbml_species.getAnnotation()
        if annotation is None:
            return []

        chebi_ids = []
        for i in range(libsbml_species.getNumCVTerms()):
            cv_term = libsbml_species.getCVTerm(i)
            if cv_term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER:
                for j in range(cv_term.getNumResources()):
                    resource = cv_term.getResourceURI(j)
                    if "chebi" in resource.lower():
                       
                        match = re.search(r'CHEBI:(\d+)', resource)
                        if match:
                            chebi_ids.append(match.group(1))

        chebi_ids = list(dict.fromkeys(chebi_ids))

        return chebi_ids
    


    @staticmethod
    def _parse_using_chebi( chebi_code: str ) -> tuple[str, int, dict]:

        """
            This function receives a string which includes ChEBI code for the compound and uses EBI API to search for the compounds chemical composition and fetches it
            Then using the chemparse package decomposes the chemical formula to its elements and returns it as a dictionary
            It should be mentioned that some ChEBI compounds don't have chemical formula registered for them and some have two different chemical compositions assigned to them

            Args:
                chebi_code (str): a 5 digit code as a string

            Returns:
                str: a string represnting the compound's chemical formula (like CH4)
                int: an integer representing the charge of the compound
                dict: A dictionary mapping molecule names (str) to their integer counts (int).
        """

        chebi_entity = chb.ChebiEntity(chebi_code)
        name = chebi_entity.get_name()

        parsed_compound={}


        formulae = chebi_entity.get_formulae() # There are two methods to fetch chemical composition of a compound: get_formulae is used when more than one composition is registerd for the compound
                                                # It returns a list containing 'Formula' classes and Getting list's length will show us that if there is one or more chemical compositions registered
        if len(formulae) == 0: # Sometimes no chemical formula is registered so the length of the list will be zero

            parsed_compound = None
            formula = None
            charge = None

        else:

            formula = chebi_entity.get_formula()
            charge = chebi_entity.get_charge()
            parsed_compound = chp.parse_formula(formula)

        return formula, charge, parsed_compound
    


    # ********************************
    # *           Function           *
    # ********************************
    def _parse_molecule_units(self, cellml_comp_code: str) -> dict:
        """
            Parses a string like '2CH4.P4' into a dictionary of molecule counts.
            For example, '2CH4.P4' becomes {'CH4': 2, 'P4': 1}.

            Args:
                cellml_comp_code (str): A string representing the encoded composition of a compound, 
                    where molecules are separated by dots ('.') and may be prefixed by a coefficient.

            Returns:
                dict: A dictionary mapping molecule names (str) to their integer counts (int).
        """

        result = defaultdict(int)
        units = cellml_comp_code.split('.')

        for unit in units:
            match = re.match(r'^(\d*)([a-zA-Z][a-zA-Z0-9]*)$', unit.strip())
            if not match:
                raise ValueError(f"Invalid format in unit: '{unit}'")
            
            qty_str, formula = match.groups()
            qty = int(qty_str) if qty_str else 1
            result[formula] += qty

        return dict(result)
    


    @staticmethod
    def _parse_compound_code(compound: str) -> tuple[int, str]:
        """
            Extracts charge and name from a compound string like '-4ATP', '2Ca', or 'H2O'.

            Args:
                compound (str): The chemical compound with its charge added to the beginning of the string

            Returns:
                str: a string which is solely the compound name
                int: the charge of the compound
        """

        match = re.match(r'^([+-]?\d*)([A-Za-z].*)$', compound)
        
        if match:
            charge_str, name = match.groups()
            
            if charge_str == '':
                charge = 0
            elif charge_str.startswith(('+', '-')):
                charge = int(charge_str)
            else:
                # number without + or - is considered positive
                charge = int(charge_str)
        else:
            # No number at the beginning → assume charge 0
            charge = 0
            name = compound

        return name, charge