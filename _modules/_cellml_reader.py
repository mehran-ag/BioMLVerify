# Importing external packages
from libcellml import Parser, Validator, Analyser, AnalyserExternalVariable, Importer, cellmlElementTypeAsString, Model as CellMLModel
import numpy as np
import chemparse as chp
import libchebipy as chb
import libsbml
import sympy as sp
import re
from collections import defaultdict

from lxml import etree

# Importing internal packages
import os
import _modules._utility as utility
from pathlib import Path, PurePath
import _modules._constants as cn

from xml.dom.minidom import parseString

from _classes.cBioMLReaction import *
from _classes.cBioMLModel import *
from _classes.cBioMLSpecies import *
from _classes.cBioMLParameter import *
from _classes.cBioMLSpeciesReference import *



class CellmlReader:



    # ********************************
    # *           Function           *
    # ********************************
    def read_file( self, file_path: Union[str, Path], cellml_strict_mode: bool = False) -> BioMLModel:
        """
            Reads a CellML file and converts it to a BioML if it is convertible.
            
            Args: 
                file_path (str or Path): The path to the CellML file.
                cellml_strict_mode (bool, optional): Whether to enforce strict CellML format rules. Defaults to False.
            
            Returns:
                BioMLModel class, to which different functions can be applied.
        """

        base_dir = os.path.dirname(file_path)

        if not isinstance(base_dir, str):
            if isinstance(base_dir, Path):
                base_dir = str(base_dir)
            elif isinstance(base_dir,PurePath):
                base_dir = str(base_dir)
            else:
                raise TypeError('The file path should be a string or a Path object.')

        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Model source file '{file_path}' does not exist.")
        
        cellml_model_name = os.path.basename(file_path)

        cellml_model = CellmlReader._read_analyse_cellml_model( file_path, cellml_strict_mode )

        cellml_contents = self._extract_cellml_content(cellml_model) #returns a dictionary containing "cellml_eqs", "cellml_flattened_eqs", "cellml_species_instances", "cellml_vars_instances", and "cellml_ast_nodes"

        cellml_vars_instances = cellml_contents["cellml_vars_instances"]

        variable_type_buckets = self._classify_variables(cellml_vars_instances)
        
        """
        "variables_type_buckets" is a dictionary as shown below. It contains all kinds of variables encoded in the CellML file

        variable_type_buckets = {
            'va': list of variables
            'co': list of coefficients
            'rc': list of rate constants
            'ra': list of rates
            'bc': list of boundary conditions
            'ev': list of equation variables
            'bv': list of boundary values
            'en': list of enzymes
        }
        """

        keys_to_check = ['va', 'co', 'rc', 'ra']

        if all(variable_type_buckets.get(k) for k in keys_to_check):

            biomlmodel_species_list = self._identify_species(variable_type_buckets['va'])

            biomlmodel_reactions_list = self._identify_reactions(variable_type_buckets['co'], biomlmodel_species_list)

            # biomlmodel_reactions_list, biomlmodel_species_list = self._identify_boundary_conditions(variable_type_buckets['bc'], biomlmodel_reactions_list, biomlmodel_species_list)

            self._find_kinetic_rate_constants(variable_type_buckets['rc'], biomlmodel_reactions_list)

            biomlmodel = BioMLModel(cellml_model.id())

            biomlmodel.species = biomlmodel_species_list

            biomlmodel.reactions = biomlmodel_reactions_list

            utility.warning_printer(f"\nThis model has been converted from a CellML model using the annotations only and not by reading the equations in the CellML file.\n")

            biomlmodel.is_direct_conversion = True

            return biomlmodel
            
        else:

            mass_action = self._find_cellml_mass_actions(**cellml_contents)

            if mass_action:

                biomlmodel = self._make_biomlmodel(cellml_model, **cellml_contents)

                return biomlmodel

            else:

                print(f"\nCellML model \"{Path(cellml_model_name).stem.upper()}\" has (an) equation(s) not governed by Mass Action Kinetics\nAnd cannot be converted to a BioML Model")

                return None






    # ********************************
    # *           Function           *
    # ********************************
    def _extract_cellml_content(self, cellml_model: CellMLModel) -> dict:
        """
            Reads a CellML file and extracts all its contents, e.g. equations and all variables.
            
            Args: 
                cellml_model: a Cellml Model which has already been read by read_file function
            
            Returns:
                dict: A dictionary containing CellML equations, variables, species, and AST nodes.
        """

        cellml_vars_instances = []

        cellml_species_instances = []

        cellml_eqs = []

        cellml_ast_nodes = []

        for i in range(cellml_model.componentCount()):

            component = cellml_model.component(i)

            num_vars = component.variableCount()

            for i in range(num_vars):

                variable = component.variable(i)

                variable_id = variable.id()

                cellml_vars_instances.append(variable)

                if variable_id:

                    cellml_species_instances.append(variable)
 

            raw_mathml = component.math()

            if raw_mathml!= '':

                if self._has_multiple_equations(raw_mathml):

                    cleaned_mathml = self._remove_units_from_mathml(raw_mathml)

                    result = self._split_equations(cleaned_mathml)

                    for i in range(len(result)):

                        ast_node = libsbml.readMathMLFromString(result[i])

                        if ast_node:

                            string_formula = libsbml.formulaToL3String(ast_node)

                            string_formula = string_formula.replace('==', '=')

                            cellml_eqs.append(string_formula)

                            cellml_ast_nodes.append(ast_node)

                else:

                    cleaned_mathml = self._remove_units_from_mathml(raw_mathml)

                    ast_node = libsbml.readMathMLFromString(cleaned_mathml)

                    if ast_node:

                        string_formula = libsbml.formulaToL3String(ast_node)

                        string_formula = string_formula.replace('==', '=')

                        cellml_eqs.append(string_formula)

                        cellml_ast_nodes.append(ast_node)

            else:

                continue

        cellml_flattened_eqs = self._flatten_equations(cellml_eqs, cellml_vars_instances)

        return {
            "cellml_eqs": cellml_eqs,
            "cellml_flattened_eqs": cellml_flattened_eqs,
            "cellml_species_instances": cellml_species_instances,
            "cellml_vars_instances": cellml_vars_instances,
            "cellml_ast_nodes": cellml_ast_nodes
        }







    # ********************************
    # *           Function           *
    # ********************************
    def _read_cellml_elements(self, cellml_model: CellMLModel, cellml_model_name: str) -> tuple[list[object], list[object]]:
        """
            Reads a CellML file and extracts all its contents, e.g. equations and all variables.
            
            Args: 
                cellml_model: a Cellml Model which has already been read by read_file function
                cellml_model_name: this is the file name (not Cellml Model's name)
            
            Returns:
                two lists:
                    list: containing CellML components (instances of cellml component class)
                    list: containing CellML variables (instances of cellml variable class)
        """

        cellml_variables = []

        cellml_components = []

        number_of_components = cellml_model.componentCount()

        if number_of_components >= 1 :

            for i in range(number_of_components):

                cellml_components.append( cellml_model.component( i ) )

        else:
            raise ValueError(f"There is no component in {cellml_model_name}")
        

        for component in cellml_components:

            number_of_variables = component.variableCount()

            for i in range(number_of_variables):

                cellml_variables.append(component.variable(i))

        return cellml_components, cellml_variables
    




    # ********************************
    # *           Function           *
    # ********************************
    def _classify_variables(self, cellml_variables: list[object]) -> dict:
        """
            Classifies the variables in a CellML Model based on their annotations:
                'va': variables
                'co': coefficients
                'rc': rate constants
                'ra': rates
                'bc': boundary conditions
                'ev': equation variables
                'bv': boundary values
                'en': enzymes
            
            Args: 
                cellml_variables: alist containing all CellML variables (instances of CellML variable class)
            
            Returns:
                dict: A dictionary containing CellML variables, coefficients, rate constants, rates, boundary conditions, equation variables, boundary values,and enzymes
        """


        variable_type_buckets = {
            'va': [],  # variables
            'co': [],  # coefficients
            'rc': [],  # rate constants
            'ra': [],  # rates
            'bc': [],  # boundary conditions
            'ev': [],  # equation variables
            'bv': [],  # boundary values
            'en': []   # enzymes
        }

        variable_classifier = {
            'va': lambda v: variable_type_buckets['va'].append(v),
            'co': lambda v: variable_type_buckets['co'].append(v),
            'rc': lambda v: variable_type_buckets['rc'].append(v),
            'ra': lambda v: variable_type_buckets['ra'].append(v),
            'bc': lambda v: variable_type_buckets['bc'].append(v),
            'ev': lambda v: variable_type_buckets['ev'].append(v),
            'bv': lambda v: variable_type_buckets['bv'].append(v),
            'en': lambda v: variable_type_buckets['en'].append(v)
        }

        for cellml_variable in cellml_variables:

            variable_id = cellml_variable.id()

            identifier = variable_id.split('_')[0]

            operation = variable_classifier.get(identifier)

            if operation:

                operation(cellml_variable)

        return variable_type_buckets





    # ********************************
    # *           Function           *
    # ********************************
    def _identify_species(self, variables: list[object]) -> list[object]:
        """
            finds all species in the model and converts them into BioML species
            
            Args: 
                cellml_variables: a list containing all CellML variables (instances of CellML variable class)
            
            Returns:
                list: a list containing all BioML species (instances of BioML species class)
        """

        biomlmodel_species_list = []

        for variable in variables:

            name = variable.name()

            matched_species = next( ( species_instance for species_instance in biomlmodel_species_list if species_instance.ID == name ), None )

            if matched_species is not None:
                continue

            cellml_id = variable.id()

            biomlmodel_species = BioMLSpecies(cellml_id)

            biomlmodel_species.name = name

            if cellml_id.split('_')[1].isdigit():

                chebi_code = cellml_id.split('_')[1]

                biomlmodel_species.chebi_code = chebi_code

                compound, charge, composition = CellmlReader._parse_using_chebi(chebi_code)

                if compound is not None:
                    biomlmodel_species.compound = compound

                if charge is not None:
                    biomlmodel_species.charge = charge

                if composition is not None:
                    biomlmodel_species.composition = composition

            else:

                compound_part = cellml_id.split('_')[1]

                if compound_part[0] == '-':
                    charge_sign = -1
                    compound_part = compound_part[1:]

                else:
                    charge_sign = +1

                compound_code = compound_part.split('-')[0]

                cellml_species_compound, cellml_species_charge = CellmlReader._parse_compound_code(compound_code)

                biomlmodel_species.charge = charge_sign * cellml_species_charge

                if len( compound_part.split('-') ) > 1:

                    cellml_species_comp_code = compound_part.split('-')[1]

                    cellml_species_composition = self._parse_molecule_units(cellml_species_comp_code)

                else:

                    cellml_species_composition = {cellml_species_compound: 1}

                biomlmodel_species.compound = cellml_species_compound

                biomlmodel_species.composition = cellml_species_composition


            biomlmodel_species_list.append(biomlmodel_species)

        return biomlmodel_species_list





    # ********************************
    # *           Function           *
    # ********************************
    def _identify_reactions(self, coefficients: list[object], biomlmodel_species_list: list[object] ) -> list[object]:
        """
            finds all reactions in the model and converts them into BioML reactions
            
            Args:
                coefficients: a list containing all coefficients (instances, which are classified as coefficients, of CellML variable class)
                biomlmodel_species_list: a list containing BioML species (intances of BioMl species class)
            
            Returns:
                list: a list containing all BioML reactions (instances of BioML reaction class)
        """

        biomlmodel_reactions_list = []

        for coefficient in coefficients:
            
            name_code = coefficient.id().split('_')[1]

            if all( char.isdigit() for char in name_code ):

                compound, _ = CellmlReader._parse_using_chebi(name_code)

            else:

                compound = name_code

            matched_biomlmodel_species =  next( ( species_instance for species_instance in biomlmodel_species_list if species_instance.compound == compound ), None )

            if matched_biomlmodel_species is None:
                raise ValueError(f"There is no match for species {compound} in the list of species")
            

            reaction_number_parts = coefficient.id().split('_')[2]

            for i, reaction_number_part in enumerate(reaction_number_parts.split('-')):

                reaction_number = reaction_number_part.split('.')[0]

                matched_biomlmodel_reaction = next( ( reaction_instance for reaction_instance in biomlmodel_reactions_list if reaction_instance.ID == reaction_number ), None )

                if matched_biomlmodel_reaction is None:
                    
                    biomlmodel_reaction = BioMLReaction(reaction_number)

                    biomlmodel_species_ref = BioMLSpeciesReference(matched_biomlmodel_species)

                    biomlmodel_species_ref.reaction_id = reaction_number

                    stoichiometry = float(coefficient.initialValue())

                    if stoichiometry > 0:

                        biomlmodel_species_ref.stoichiometry = abs(stoichiometry)

                        biomlmodel_reaction.products.append(biomlmodel_species_ref)

                    elif stoichiometry < 0:

                        biomlmodel_species_ref.stoichiometry = abs(stoichiometry)

                        biomlmodel_reaction.reactants.append(biomlmodel_species_ref)

                    else:

                        raise ValueError(f"The value of stoichiometric coefficient is not acceptable: {coefficient.initialValue()}")

                    biomlmodel_reactions_list.append(biomlmodel_reaction)

                else:

                    if i > 0:
                        continue

                    biomlmodel_species_ref = BioMLSpeciesReference(matched_biomlmodel_species)

                    biomlmodel_species_ref.reaction_id = reaction_number

                    stoichiometry = float(coefficient.initialValue())

                    if stoichiometry > 0:

                        biomlmodel_species_ref.stoichiometry = abs(stoichiometry)

                        matched_biomlmodel_reaction.products.append(biomlmodel_species_ref)

                    elif stoichiometry < 0:

                        biomlmodel_species_ref.stoichiometry = abs(stoichiometry)

                        matched_biomlmodel_reaction.reactants.append(biomlmodel_species_ref)

                    else:

                        raise ValueError(f"The value of stoichiometric coefficient is not acceptable: {coefficient.initialValue()}")
                    
        return biomlmodel_reactions_list







    # ********************************
    # *           Function           *
    # ********************************
    def _identify_boundary_conditions(self, boundary_conditions: list[object], biomlmodel_reactions_list: list[object], biomlmodel_species_list: list[object]) -> tuple[list[object], list[object]]:
        """
            gets boundary conditions in the model and converts them into equivalent BioML reactions
            
            Args:
                boundary_conditions: a list containing all variables annotated as boundary conditions
                biomlmodel_reactions_list: a list containing BioML reactions (intances of BioMl reaction class)
                biomlmodel_species_list: a list containing BioML species (intances of BioMl species class)
            
            Returns:
                two lists:
                    list: containing BioML reactions (instances of BioML reaction class)
                    list: containing BioML species (instances of BioML species class)
        """

        if boundary_conditions:

            for bc in boundary_conditions:

                reaction_number = bc.id()

                if not any(biomlmodel_reaction.ID == reaction_number for biomlmodel_reaction in biomlmodel_reactions_list):
                        
                    biomlmodel_reaction = BioMLReaction(reaction_number)

                    biomlmodel_reaction.boundary_condition = True

                    name_code =  bc.id().split('_')[1].split('.')[0]

                    try:

                        flow_direction = bc.id().split('_')[2].split('.')[0]

                    except IndexError:

                        flow_direction = ''

                    if all( char.isdigit() for char in name_code ):

                        compound, _ = CellmlReader._parse_using_chebi(name_code)

                    else:

                        compound = name_code


                    matched_species = next( ( species_instance for species_instance in biomlmodel_species_list if species_instance.compound == compound ), None )

                    if matched_species is None:
                        raise ValueError(f"The Species for the boundary condition {reaction_number} can't be found!")
                    
                    
                    biomlmodel_species_ref = BioMLSpeciesReference(matched_species)

                    biomlmodel_species_ref.reaction_id = reaction_number

                    biomlmodel_species_ref.stoichiometry = 1

                    if flow_direction == 'i':

                        biomlmodel_reaction.products.append(biomlmodel_species_ref)

                    elif flow_direction == 'o':

                        biomlmodel_reaction.reactants.append(biomlmodel_species_ref)

                    else:

                        biomlmodel_reaction.reactants.append(biomlmodel_species_ref)
                    

                    # Now, we will create the virtual internal and external species for the compound
                    compound_bc = compound + '_e'

                    matched_bc_species =  next( ( species_instance for species_instance in self._biomlmodel_species_list if species_instance.compound == compound_bc ), None )

                    if matched_bc_species is None:

                        biomlmodel_species = BioMLSpecies(compound_bc)

                        biomlmodel_species.compound = compound_bc

                        biomlmodel_species_list.append(biomlmodel_species)

                    else:

                        biomlmodel_species = matched_bc_species

                    biomlmodel_species_ref = BioMLSpeciesReference(biomlmodel_species)

                    biomlmodel_species_ref.reaction_id = reaction_number

                    biomlmodel_species_ref.stoichiometry = 1

                    if flow_direction == 'i':

                        biomlmodel_reaction.reactants.append(biomlmodel_species_ref)

                    elif flow_direction == 'o':

                        biomlmodel_reaction.products.append(biomlmodel_species_ref)

                    else:

                        biomlmodel_reaction.products.append(biomlmodel_species_ref)

                    biomlmodel_reactions_list.append(biomlmodel_reaction)

        return biomlmodel_reactions_list, biomlmodel_species_list






    # ********************************
    # *           Function           *
    # ********************************
    def _find_kinetic_rate_constants(self, rate_constants: list[object], biomlmodel_reactions_list: list[object]) -> None:
        """
            reads equations stored for a reaction and indetifies reaction rate constants for the reaction and updates the BioML reaction instance
            
            Args:
                rate_constants: a list containing all variables annotated as reaction rate constants
                biomlmodel_reactions_list: a list containing BioML reactions (intances of BioMl reaction class)
            
            Returns:
               None
        """

        for rate_constant in rate_constants:

            rc_id = rate_constant.id()

            reaction_number_parts = rc_id.split('_')[2]

            read_reactions = []

            for reaction_number_part in reaction_number_parts.split('-'):

                reaction_name = reaction_number_part.split('.')[0]

                if reaction_name in read_reactions:
                    continue
                else:
                    read_reactions.append(reaction_name)

                direction = rc_id.split('_')[1]

                matched_reaction =  next( ( reaction_instance for reaction_instance in biomlmodel_reactions_list if reaction_instance.ID == reaction_name ), None )

                if matched_reaction is None:
                    raise ValueError(f"Reaction {reaction_name} cannot be found in the reaction list")
                
                if direction == "f":
                    matched_reaction.kinetic_forward_rate_constant = rate_constant.name()
                    matched_reaction.kinetic_forward_rate_constant_value = float(rate_constant.initialValue())

                elif direction == "r":
                    matched_reaction.kinetic_reverse_rate_constant = rate_constant.name()
                    matched_reaction.kinetic_reverse_rate_constant_value = float(rate_constant.initialValue())

                    matched_reaction.reversible = True





    # ********************************
    # *           Function           *
    # ********************************
    def _find_cellml_mass_actions(self, **cellml_contents) -> bool:
        """
            checks the equations for the patterns similar to Mass Action formulations and marks the equation as Mass Action if pattern is found in the reaction
            
            Keyword Args:
                cellml_vars_instances: a list containing CellML variables (instances of variable class)
                cellml_ast_nodes: a list containing Abstract Syntax Trees (AST), a definition of libsbml specific class
                cellml_flattened_eqs : a list containing all flattened equations in the CellML model (All variables are replaced by equations if they are representative of any equation in Flattened equations)
                species_in_cellml_eqs: a list containing all CellML variables annotated as species used in CellML equations
            
            Returns:
               bool: True if any mass action pattern is found, False otherwise.
        """

        cellml_vars_instances = cellml_contents["cellml_vars_instances"]

        cellml_ast_nodes = cellml_contents["cellml_ast_nodes"]

        cellml_flattened_eqs = cellml_contents["cellml_flattened_eqs"]

        species_in_cellml_eqs = self._find_species_in_all_eqs(**cellml_contents)

        mass_action = True

        cellml_vars = [var_instance.name() for var_instance in cellml_vars_instances]

        for cellml_ast_node in cellml_ast_nodes:

            if cellml_ast_node:

                cellml_eq = libsbml.formulaToL3String(cellml_ast_node)

                if '==' in cellml_eq:

                    cellml_eq_lhs = cellml_eq.split('==')[0].strip()

                if cellml_flattened_eqs.get(cellml_eq_lhs) is not None:
                    cellml_eq_rhs = str(cellml_flattened_eqs[cellml_eq_lhs])

                    simp_cellml_eq = str(sp.simplify(cellml_flattened_eqs[cellml_eq_lhs]))

                    eq_mass_action = self._check_mass_action( cellml_eq_rhs, simp_cellml_eq, cellml_vars, species_in_cellml_eqs )

                    if not eq_mass_action:

                        mass_action = eq_mass_action

        return mass_action


    

    # ********************************
    # *           Function           *
    # ********************************
    def _find_species_in_all_eqs(self, **cellml_contents) -> list[str]:
        """
            checks the equations for the existence of species and stores the species in a list
            
            Keyword Args:
                cellml_species_instances: a list containing CellML variables annotated as species (instances of variable class)
                cellml_ast_nodes: a list containing Abstract Syntax Trees (AST), a definition of libsbml specific class
            
            Returns:
               list: containing CellMl variables used as species in all equations of the CellML model
        """

        cellml_ast_nodes = cellml_contents["cellml_ast_nodes"]

        cellml_species_instances = cellml_contents["cellml_species_instances"]

        cellml_species = [species_instance.name() for species_instance in cellml_species_instances]

        species_in_cellml_eqs = []

        for cellml_ast_node in cellml_ast_nodes:

            if cellml_ast_node:

                cellml_eq_vars = CellmlReader._get_variables(cellml_ast_node)

                for cellml_eq_var in cellml_eq_vars:

                    if cellml_eq_var in cellml_species:

                        species_in_cellml_eqs.append(cellml_eq_var)

        species_in_cellml_eqs = list(dict.fromkeys(species_in_cellml_eqs))

        return species_in_cellml_eqs





    # ********************************
    # *           Function           *
    # ********************************
    def _check_mass_action(self, cellml_eq: str, simp_cellml_eq: str, cellml_vars: list[str], species_in_cellml_eq: list[str]) -> bool:
        """
            checks a single equation for the patterns found in Mass Action Kinetics equations
            
            Keyword Args:
                cellml_eq (str): a CellML equations as a string
                simp_cellml_eq (str): a simplified (using sympy simplification command) CellML equation as a string
                cellml_vars (list): a list containing all CellML variables
                species_in_cellml_eq (list): a list containing the species used in this specific equation
            
            Returns:
               bool: True if any mass action pattern is found, False otherwise.
        """

        flag = False

        if self._single_product( cellml_eq, simp_cellml_eq ) or \
            self._diff_of_products( cellml_eq, simp_cellml_eq ):

                if len( species_in_cellml_eq ) != 0:

                    flag = True

                    if len( species_in_cellml_eq ) == 1:

                        if cellml_eq.count(species_in_cellml_eq[0]) != 1 or simp_cellml_eq.count(species_in_cellml_eq[0]) != 1:
                            flag = False
        
        try:
            fracs = self._find_frac_parts(simp_cellml_eq, cellml_vars)
            if len(species_in_cellml_eq) > 0:
                for i in range(len(species_in_cellml_eq)):
                    if species_in_cellml_eq[i] in fracs["denominator"]:
                        flag = False

        except:
            pass

        return flag




    # ********************************
    # *           Function           *
    # ********************************
    def _single_product(self, kinetic_law: str, simple_kinetic_law: str) -> bool:
        """
            checks the input kinetic law to see if it is a single product of terms
            
            Keyword Args:
                kinetic_law (str): kinetic law of a reaction
                simple_kinetic_law (str): a simplified (using sympy simplification command) kinetic law of a reaction
            
            Returns:
               bool: True if is single product of terms, False otherwise.
        """

        flag = True

        if "+" in kinetic_law or "-" in kinetic_law:
            flag = False
            if "e-" in kinetic_law or "exp(-" in kinetic_law:
                flag = True
        elif "+" in simple_kinetic_law or "-" in simple_kinetic_law:
            flag = False
            if "e-" in simple_kinetic_law or "exp(-" in simple_kinetic_law:
                flag = True

        return flag


    # ********************************
    # *           Function           *
    # ********************************
    def _diff_of_products(self, kinetic_law: str, simple_kinetic_law: str) -> bool:
        """
            checks the input kinetic law to see if it is difference of product of two terms
            
            Keyword Args:
                kinetic_law (str): kinetic law of a reaction
                simple_kinetic_law (str): a simplified (using sympy simplification command) kinetic law of a reaction
            
            Returns:
               bool: True if is difference of product of two terms, False otherwise.
        """

        flag = False

        if self._single_product(kinetic_law, simple_kinetic_law) == False:
            terms = kinetic_law.split("-")
            if len(terms) == 2:
                flag = True

        return flag
    

    # ********************************
    # *           Function           *
    # ********************************
    def _find_frac_parts( self, simp_cellml_eq: str, cellml_vars: list[str] ) -> dict:
        """
            Finds the numerator and denominator of a fraction in a kinetic law string.

            Keyword Args:
                simp_cellml_eq (str): Simplified CellMl equation using Sympy simplification method
                cellml_vars (list[str]): A list containing CellML variables as strings

            Returns:
                dict: A dictionary with the following keys:
                    - "numerator": A string containing the numerator of the fraction.
                    - "denominator": A string containing the denominator of the fraction.
        """

        fracs = {"numerator": '', "denominator": ''}

        try:
            symbol_dict = {cellml_var: sp.symbols(cellml_var) for cellml_var in cellml_vars}
            cellml_sympy_eq = sp.sympify(simp_cellml_eq, locals=symbol_dict)
            numerator, denominator = cellml_sympy_eq.as_numer_denom()
            fracs["numerator"] = str(numerator)
            fracs["denominator"] = str(denominator)
        except Exception:
            pass

        return fracs
    


    


    # ********************************
    # *           Function           *
    # ********************************
    def _has_multiple_equations(self, mathml_str: str) -> int:
        """
            Checks whether there are more than one <apply><eq/>...</apply> elements in the MathML.
            Each <apply><eq/>...</apply> means an equation in MathML, so this function finds the number of equations stored in a mathml script.

            Args:
                mathml_str (str): a string of a mathml script

            Returns:
                int: the number of equations in the mathml script
        """

        root = etree.fromstring(mathml_str.encode())

        NSMAP = {'m': 'http://www.w3.org/1998/Math/MathML'}

        eq_applies = root.xpath(".//m:apply[m:eq]", namespaces=NSMAP)

        return len(eq_applies) > 1
    





    # ********************************
    # *           Function           *
    # ********************************
    def _split_equations(self, mathml_str: str) -> list:
        """
            Splits a MathML script containing multiple equations into individual <math> blocks.

            Args:
                mathml_str (str): a string of mathml script

            Returns:
                list: a list of MathML strings representing an equation
        """
        root = etree.fromstring(mathml_str.encode())
        NSMAP = {'m': 'http://www.w3.org/1998/Math/MathML'}
        equations = root.xpath(".//m:apply[m:eq]", namespaces=NSMAP)

        result = []
        for eq in equations:
            # Create a new <math> element
            MATHML_NS = "http://www.w3.org/1998/Math/MathML"
            math_elem = etree.Element("math")
            math_elem.set("xmlns", MATHML_NS)
            # Append a deep copy of the equation node
            math_elem.append(eq)
            mathml = etree.tostring(math_elem, pretty_print=True, xml_declaration=False).decode()
            result.append(mathml.strip())
        return result
    





    # ********************************
    # *           Function           *
    # ********************************
    def _remove_units_from_mathml(self, mathml_str: str) -> str:
        """
            Finds attributes for a number element in mathml script and removes it.
            The existence of the unit confuses libsbml function converting mathml script into a string

            Args:
                mathml_str (str): a string of mathml script

            Returns:
                str: a string of mathml script having its units removed
        """

        root = etree.fromstring(mathml_str.encode())

        ns = {
            "m": "http://www.w3.org/1998/Math/MathML",
        }

        for cn in root.xpath(".//m:cn", namespaces=ns):
            cn.attrib.pop('{http://www.cellml.org/cellml/1.0#}units', None)
            cn.attrib.pop('{http://www.cellml.org/cellml/2.0#}units', None)

        return etree.tostring(root).decode()




    # ********************************
    # *           Function           *
    # ********************************
    def _flatten_equations(self, cellml_eqs: list[str], cellml_vars_instances: list[object]) -> list:
        """
            Returns a dictionary mapping variables (as strings) to flattened equations (sympy expressions) where all variables defined by equations have been substituted by their defnitions.

            Args:
                cellml_eqs (list): containing all equations as they are imported from CellML
                cellml_vars_instances (list): containing all variables in a CellML model

            Returns:
                list: containing the flattened equations
        """

        # Step 1: Create symbol dictionary
        cellml_vars = [var_instance.name() for var_instance in cellml_vars_instances]
        symbol_dict = {var: sp.symbols(var) for var in cellml_vars}

        # Step 2: Build equation dictionary: LHS string -> RHS sympy expression
        eq_dict = {}
        for cellml_eq in cellml_eqs:
            if '=' in cellml_eq:
                lhs_str = cellml_eq.split('=')[0].strip()
                rhs_str = cellml_eq.split('=')[1].strip()
                rhs_expr = sp.sympify(rhs_str.replace("^", "**"), locals=symbol_dict)
                eq_dict[lhs_str] = rhs_expr

        # Step 3: Recursive substitution
        def substitute_all(expr, eq_dict):
            prev_expr = None
            while expr != prev_expr:
                prev_expr = expr
                for var_str, sub_expr in eq_dict.items():
                    expr = expr.subs(symbol_dict[var_str], sub_expr)
            return expr

        # Step 4: Determine which variables are used as intermediate
        substituted_vars = set()
        for rhs_expr in eq_dict.values():
            for symbol in rhs_expr.free_symbols:
                var_name = str(symbol)
                if var_name in eq_dict:
                    substituted_vars.add(var_name)

        # Step 5: Flatten and filter
        flattened_eqs = {}
        for var, rhs_expr in eq_dict.items():
            if var not in substituted_vars:
                flattened_rhs = substitute_all(rhs_expr, eq_dict)
                flattened_eqs[var] = sp.simplify(flattened_rhs)

        return flattened_eqs
    




    # ********************************
    # *           Function           *
    # ********************************
    def _make_biomlmodel(self, cellml_model: CellMLModel, **cellml_contents) -> BioMLModel:

        """
            Returns a dictionary mapping variables (as strings) to flattened equations (sympy expressions) where all variables defined by equations have been substituted by their defnitions.

            Args:
                cellml_model (CellMLModel): CellML model
                cellml_contents( dict): a dictionary containing {"cellml_eqs": cellml_eqs,
                    "cellml_flattened_eqs": cellml_flattened_eqs,
                    "cellml_species_instances": cellml_species_instances,
                    "cellml_vars_instances": cellml_vars_instances,
                    "cellml_ast_nodes": cellml_ast_nodes}

            Returns:
               BioMLModel
        """


        species_in_cell_equations = self._find_species_in_all_eqs(**cellml_contents)

        cellml_flattened_eqs = cellml_contents["cellml_flattened_eqs"]

        cellml_species_instances = cellml_contents["cellml_species_instances"]

        cellml_vars_instances = cellml_contents["cellml_vars_instances"]

        reversible = True

        if cellml_model.id() != '':

            biomlmodel_name = cellml_model.id()

        elif cellml_model.name() != '':

            biomlmodel_name = cellml_model.name()

        else:

            biomlmodel_name = 'Not assigned'

        biomlmodel = BioMLModel(biomlmodel_name)

        biomlmodel_reactions_list = []

        for reaction_id, sympy_expr_eq in cellml_flattened_eqs.items():

            rates = self._get_forward_reverse_rate_expressions(sympy_expr_eq)

            if len(rates["forward_rate"]) == 0:
                raise ValueError(f"No forward reaction rate has been found in equation: {str(sympy_expr_eq)}")

            elif len(rates["forward_rate"]) > 1:
                raise ValueError(f"There are more than one addition terms in the forward reaction rate of the kinetic law equation: {str(sympy_expr_eq)}")
            
            else:

                forward_rate_expr = rates["forward_rate"][0]

                forward_rate_contents = self._analyze_sympy_expression(forward_rate_expr, species_in_cell_equations)

            if len(rates["reverse_rate"]) == 0:
                reversible = False

            elif len(rates["reverse_rate"]) == 1:

                reverse_rate_expr = rates["reverse_rate"][0]

                reverse_rate_contents = self._analyze_sympy_expression(reverse_rate_expr, species_in_cell_equations)

            else:

                raise ValueError(f"There are more than one addition terms in the reverse reaction rate of the kinetic law equation: {str(sympy_expr_eq)}")
            
            if not reversible:
                raise ValueError(f"Reaction {reaction_id} is not reversible and the reaction can not be defined from the reaction rate equation")


            biomlmodel_reaction = BioMLReaction(reaction_id)

            biomlmodel_reaction.kinetic_law = str(sympy_expr_eq).replace('**', '^')

            biomlmodel_reaction.kinetic_law_type = 'Mass Action'

            biomlmodel_species_list = []

            biomlmodel_reactants_list = []

            biomlmodel_products_list = []

            local_parameters = []

            forward_rate_constant = str(forward_rate_contents["rate_constant"])

            cellml_forward_rate_instance = self._return_cellml_parameter_instance( forward_rate_constant, cellml_vars_instances )

            if cellml_forward_rate_instance is not None:

                parameter_name = cellml_forward_rate_instance.name()

                parameter_value = float(cellml_forward_rate_instance.initialValue())

                biomlmodel_reaction.kinetic_forward_rate_constant = forward_rate_constant

                biomlmodel_reaction.kinetic_forward_rate_constant_value = parameter_value

                biomlmodel_parameter = BioMLParameter(parameter_name)

                biomlmodel_parameter.value = parameter_value

                local_parameters.append(biomlmodel_parameter)

            for species_name, stoichiometry in forward_rate_contents["stoichiometry"].items():

                cellml_species_instance = self._return_cellml_species_instance( str(species_name), cellml_species_instances )

                if cellml_species_instance is not None:

                    cellml_id = cellml_species_instance.id()

                    biomlmodel_species = BioMLSpecies(cellml_id)

                    biomlmodel_species.name = str(species_name)

                    if cellml_id.split('_')[1].isdigit():

                        chebi_code = cellml_id.split('_')[1]

                        biomlmodel_species.chebi_code = chebi_code

                        compound, charge, composition = CellmlReader._parse_using_chebi(chebi_code)

                        if compound is not None:
                            biomlmodel_species.compound = compound

                        if charge:
                            biomlmodel_species.charge = charge

                        if composition is not None:
                            biomlmodel_species.composition = composition

                    else:

                        compound_part = cellml_id.split('_')[1]

                        if compound_part[0] == '-':
                            charge_sign = -1
                            compound_part = compound_part[1:]

                        else:
                            charge_sign = +1

                        compound_code = compound_part.split('-')[0]

                        cellml_species_compound, cellml_species_charge = CellmlReader._parse_compound_code(compound_code)

                        biomlmodel_species.charge = charge_sign * cellml_species_charge

                        if len( compound_part.split('-') ) > 1:

                            cellml_species_comp_code = compound_part.split('-')[1]

                            cellml_species_composition = self._parse_molecule_units(cellml_species_comp_code)

                        else:

                            cellml_species_composition = {cellml_species_compound: 1}
                            

                        biomlmodel_species.compound = cellml_species_compound

                        biomlmodel_species.composition = cellml_species_composition


                    biomlmodel_species_list.append(biomlmodel_species)

                    biomlmodel_species_reference = BioMLSpeciesReference(biomlmodel_species)

                    biomlmodel_species_reference.reaction_id = reaction_id

                    biomlmodel_species_reference.stoichiometry = stoichiometry

                    biomlmodel_reactants_list.append(biomlmodel_species_reference)

            reverse_rate_constant = str(reverse_rate_contents["rate_constant"])

            cellml_reverse_rate_instance = self._return_cellml_parameter_instance(reverse_rate_constant, cellml_vars_instances)

            parameter_id = cellml_reverse_rate_instance.id()

            parameter_value = float(cellml_reverse_rate_instance.initialValue())

            biomlmodel_reaction.kinetic_reverse_rate_constant = reverse_rate_constant

            biomlmodel_reaction.kinetic_reverse_rate_constant_value = parameter_value

            biomlmodel_parameter = BioMLParameter(parameter_id)

            biomlmodel_parameter.value = parameter_value

            local_parameters.append(biomlmodel_parameter)

            for species_name, stoichiometry in reverse_rate_contents["stoichiometry"].items():

                cellml_species_instance = self._return_cellml_species_instance( str(species_name), cellml_species_instances )

                if cellml_species_instance is not None:

                    cellml_id = cellml_species_instance.id()

                    biomlmodel_species = BioMLSpecies(cellml_id)

                    biomlmodel_species.name = str(species_name)

                    if cellml_id.split('_')[1].isdigit():

                        chebi_code = cellml_id.split('_')[1]

                        biomlmodel_species.chebi_code = chebi_code

                        compound, charge, composition = CellmlReader._parse_using_chebi(chebi_code)

                        if compound is not None:
                            biomlmodel_species.compound = compound

                        if charge:
                            biomlmodel_species.charge = charge

                        if composition is not None:
                            biomlmodel_species.composition = composition

                    else:

                        compound_part = cellml_id.split('_')[1]

                        if compound_part[0] == '-':
                            charge_sign = -1
                            compound_part = compound_part[1:]

                        else:
                            charge_sign = +1

                        compound_code = compound_part.split('-')[0]

                        cellml_species_compound, cellml_species_charge = CellmlReader._parse_compound_code(compound_code)

                        biomlmodel_species.charge = charge_sign * cellml_species_charge

                        if len( compound_part.split('-') ) > 1:

                            cellml_species_comp_code = compound_part.split('-')[1]

                            cellml_species_composition = self._parse_molecule_units(cellml_species_comp_code)

                        else:

                            cellml_species_composition = {cellml_species_compound: 1}

                        biomlmodel_species.compound = cellml_species_compound

                        biomlmodel_species.composition = cellml_species_composition

                    biomlmodel_species_list.append(biomlmodel_species)

                    biomlmodel_species_reference = BioMLSpeciesReference(biomlmodel_species)

                    biomlmodel_species_reference.reaction_id = reaction_id

                    biomlmodel_species_reference.stoichiometry = stoichiometry

                    biomlmodel_products_list.append(biomlmodel_species_reference)

            for species in biomlmodel_species_list:
                if species not in biomlmodel.species:
                    biomlmodel.species.append(species)

            biomlmodel_reaction.reversible = reversible

            biomlmodel_reaction.local_parameters = local_parameters

            biomlmodel_reaction.reactants = biomlmodel_reactants_list

            biomlmodel_reaction.products = biomlmodel_products_list

            biomlmodel_reactions_list.append(biomlmodel_reaction)

        biomlmodel.reactions = biomlmodel_reactions_list

        biomlmodel.species = biomlmodel_species_list

        return biomlmodel

        


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


    

    def _return_cellml_species_instance(self, species_name: str, cellml_species_instances: list) -> object:
        """
            Returns the CellML variable instance whose name matches the given species_name.

            Args:
                species_name (str): The species name as a string.
                cellml_species_instances (list): A list of CellML species instances (instances of the CellML variable class).

            Returns:
                object: The matching CellML variable instance.
        """

        matched_instance = None

        for cellml_species_instance in cellml_species_instances:

            if cellml_species_instance.name() == species_name:

                matched_instance = cellml_species_instance


        return matched_instance
    

    def _return_cellml_parameter_instance(self, var_name: str, cellml_vars_instances: list) -> object:
        """
            Returns the CellML variable instance whose name matches the given var_name.

            Args:
                var_name (str): The parameter name as a string.
                cellml_vars_instances (list): A list of CellML variable instances (instances of the CellML variable class).

            Returns:
                object: The matching CellML variable instance.
        """

        matched_instance = None

        for cellml_var_instance in cellml_vars_instances:

            if cellml_var_instance.name() == var_name:

                matched_instance = cellml_var_instance


        return matched_instance



    # ********************************
    # *           Function           *
    # ********************************
    def _process_rate_expression(self, rate_expr: str, species: list) -> dict:
        """
            Parses a rate expression to identify species with their stoichiometric powers,
            and modifies the expression by replacing species terms with '1'.

            Args:
                rate_expr (str): The original rate expression as a string.
                species (list): A list of species names (strings) to look for in the expression.

            Returns:
                dict: A dictionary containing:
                    - 'modified_expression' (str): The rate expression after replacing species terms with '1'.
                    - 'stoichiometry' (dict): A mapping from species names to their stoichiometric powers (integers).
        """

        stoichiometry = {}
        expr = rate_expr  # Make a copy to modify

        for sp in sorted(species, key=len, reverse=True):
            # First, handle exponentiated species (e.g., x_NO**2)
            power_match = re.search(rf"{re.escape(sp)}\s*\*\*\s*(\d+)", expr)
            if power_match:
                stoichiometry[sp] = int(power_match.group(1))
                expr = re.sub(rf"{re.escape(sp)}\s*\*\*\s*\d+", "1", expr)
            else:
                # Then handle species without exponent
                pattern_no_power = rf"(?<![\w]){re.escape(sp)}(?![\w])"
                match = re.search(pattern_no_power, expr)
                if match:
                    stoichiometry[sp] = 1
                    expr = re.sub(pattern_no_power, "1", expr)

        return {
            'modified_expression': expr,
            'stoichiometry': stoichiometry
        }
    




    # ********************************
    # *           Function           *
    # ********************************
    def _analyze_sympy_expression(self, rate_expr: sp.Expr, species_str_list: list[str]) -> dict:
        """
            Analyzes a SymPy rate expression to extract stoichiometric powers of species
            and isolate the constant rate component.

            Args:
                rate_expr (sympy.Basic): A SymPy expression representing the rate law.
                species_str_list (list): A list of species names (as strings) present in the expression.

            Returns:
                dict: A dictionary containing:
                    - 'stoichiometry' (dict): A mapping of species names to their integer exponents.
                    - 'rate_constant' (sympy.Basic): The remaining expression after substituting species with 1.

            Raises:
                ValueError: If the expression is not a product of terms (i.e., not a `Mul` expression).
        """
        
        # Convert species names (strings) to sympy symbols
        species_syms = [sp.Symbol(name) for name in species_str_list]

        # Get species and their exponents from the expression
        if rate_expr.is_Mul:
            powers = rate_expr.as_powers_dict()
        else:
            raise ValueError("Rate expression is not product of terms, so powers of variables cannot be extracted by Sympy!")
        
        stoichiometry = {}

        for sp_sym, sp_name in zip(species_syms, species_str_list):
            power = powers.get(sp_sym)
            if power is not None:
                stoichiometry[sp_name] = int(power)

        # Replace species with 1 to isolate the constant part
        substituted_expr = rate_expr.subs({sym: 1 for sym in species_syms})

        return {
            'stoichiometry': stoichiometry,
            'rate_constant': substituted_expr
        }






    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # *      Internal Function       *
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    def _get_forward_reverse_rate_expressions(self, expression: sp.Expr) -> dict:
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
                    rates["reverse_rate"].append(-argument)
                else:
                    rates["forward_rate"].append(argument)

        else:
            rates["forward_rate"].append(expression)

        return rates
    # --------------------------------------------------
    # --------------------------------------------------






    # ********************************
    # *           Function           *
    # ********************************
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

            charge = chebi_entity.get_charge()
            formula = chebi_entity.get_formula()
            parsed_compound = chp.parse_formula(formula)

        return formula, charge, parsed_compound






    # ********************************
    # *           Function           *
    # ********************************
    @staticmethod
    def _get_formula_using_chebi ( chebi_code: str ) -> str:
    
        """
            This function receives a string which includes ChEBI code for the compound and uses EBI API to search for the compound's chemical composition and fetches it
            It should be mentioned that some ChEBI compounds don't have chemical formula registered for them and some have two different chemical compositions assigned to them

            Args:
                chebi_code (str): a 5 digit code as a string

            Returns:
                str: a string represnting the compound's chemical formula (like CH4)
        """

        chebi_entity = chb.ChebiEntity(chebi_code)

        formulae = chebi_entity.get_formulae() # There are two methods to fetch chemical composition of a compound: get_formulae is used when more than one composition is registerd for the compound
                                                # It returns a list containing 'Formula' classes and Getting list's length will show us that if there is one or more chemical compositions registered
        if len(formulae) == 0: # Sometimes no chemical formula is registered so the length of the list will be zero

            formula = None

        elif len(formulae) == 1:

            formula = chebi_entity.get_formula()


        else:

            formula = formulae[1].get_formula()


        return formula
    

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





    # ********************************
    # *           Function           *
    # ********************************
    @staticmethod
    def _read_analyse_cellml_model( file_path: str, cellml_strict_mode: bool ) -> CellMLModel:
        """
            Reads, validates, resolves, and analyses a CellML model from a given file.

            This function performs the following steps on a CellML model file:
            1. Parses the model using the libCellML parser.
            2. Validates the model structure and syntax.
            3. Resolves imports and flattens the model into a single, self-contained structure.
            4. Adds external variables and their dependencies.
            5. Analyses the model for any semantic or structural issues.

            If any critical issues are found during parsing, validation, import resolution,
            or flattening, an appropriate `ValueError` is raised. Non-critical warnings
            are printed.

            Parameters:
                file_path (str): The full path to the CellML file to be read and analysed.
                cellml_strict_mode (bool): Whether to parse and import the model using strict mode.
                                        When True, the parser enforces stricter compliance rules.

            Returns:
                CellMLModel: The successfully parsed, validated, and analysed CellML model.

            Raises:
                TypeError: If the file cannot be parsed properly due to syntax issues.
                ValueError: If the model fails validation, import resolution, or flattening,
                            or if external variables cannot be added.

            Notes:
                - All warnings during parsing, validation, and import resolution are printed
                to the console using a `utility.message_printer`.
                - The function uses the libCellML `Parser`, `Validator`, `Importer`, and `Analyser`
                interfaces to process the model.
                - External variables and their dependencies are managed via an analyser
                before the model is analysed.
        """

        model_name = os.path.basename(file_path)
            
        parser = Parser(cellml_strict_mode)

        with open(file_path, 'r') as f:
            cellml_model = parser.parseModel(f.read())

            no_parser_warnings = parser.issueCount()

            if no_parser_warnings > 0:

                for i in range(no_parser_warnings):
                    utility.message_printer(parser.issue(i).description())

        validator = Validator()
        validator.validateModel(cellml_model)

        no_validator_warnings = validator.issueCount()

        if no_validator_warnings > 0:

            for i in range(no_validator_warnings):
                utility.message_printer(validator.issue(i).description(), color="yellow")

            raise ValueError(f"Model {model_name} has validation issues and cannot be imported!")


        importer = Importer(cellml_strict_mode)

        base_dir = os.path.dirname(file_path)
        importer.resolveImports(cellml_model, base_dir)

        no_importer_warnings = importer.issueCount()

        if no_importer_warnings > 0:

            for i in range(no_importer_warnings):
                utility.message_printer( importer.issue(i).description(), color="yellow")

            raise ValueError(f"Model {model_name} has import issues and cannot be imported!")

        else:
            if cellml_model.hasUnresolvedImports():
                utility.message_printer("There are Unresolved Import Issues", color="magenta")

                raise ValueError(f"Model {model_name} has import issues and cannot be imported!")


        flat_cellml_model = importer.flattenModel(cellml_model)

        if not flat_cellml_model:
            raise ValueError(f"Model {model_name} cannot be imported: flattening issues!")
        
        external_variables_dic = CellmlReader._ext_var_dic(flat_cellml_model)

        analyser = Analyser()

        for external_variable in external_variables_dic.keys():
            aev=AnalyserExternalVariable(external_variable)
            for dependency in external_variables_dic[external_variable]:
                aev.addDependency(dependency)
            if not analyser.addExternalVariable(aev):
                raise ValueError(f"External variables cannot be added to the model in {model_name}")
            
        analyser.analyseModel(flat_cellml_model)

        no_analyser_warnings = analyser.issueCount()

        if no_analyser_warnings > 0:

            for i in range(no_analyser_warnings):
                utility.message_printer("\n" + analyser.issue(i).description(), color="yellow")
        

        return cellml_model








    # ********************************
    # *           Function           *
    # ********************************
    @staticmethod
    def _ext_var_dic(flatModel: CellMLModel, external_variables_info: dict = None):
        """
            Create a dictionary of external variables in the flattened model.
            
            Args
                flatModel: Model
                    The flattened CellML model.
                external_variables_info: dict
                    The external variables to be specified, in the format of {id:{'component': , 'name': }}
            
            Returns:
                dict: The dictionary of external variables in the flattened model, in the format of {external_variable:[]}

            Raises:
                ValueError: If an external variable is not found in the flattened model.
            
            Notes:
                No dependency is specified for the external variables.
        """

        if external_variables_info is None:
            external_variables_info = {}

        external_variables_dic={}
        for _, ext_var_info in external_variables_info.items():
            flat_ext_var=flatModel.component(ext_var_info['component']).variable(ext_var_info['name']) # TODO: check if equivalent variables should be considered
            if flat_ext_var:
                external_variables_dic[flat_ext_var]=[]
            else:
                raise ValueError ("The external variable {} in component {} is not found in the flattened model!"
                                .format(ext_var_info['component'],ext_var_info['name']))

        return external_variables_dic
    


    @staticmethod
    def _identify_variables(ast_node: libsbml.ASTNode, result: list[str]) -> list[str]:
        """
            Recursively traverses an SBML AST (Abstract Syntax Tree) to identify and collect variable names.

            This method inspects the children of a given `ASTNode` and recursively extracts variable names,
            skipping function nodes and continuing traversal for unnamed nodes or functions. It accumulates
            results into the provided `result` list and returns it.

            Parameters:
                ast_node (libsbml.ASTNode): The current AST node to process.
                result (list): A list to accumulate identified variable names during recursion.

            Returns:
                list: A list of variable names (strings) extracted from the AST.

            Raises:
                ValueError: If the recursion depth exceeds `cn.MAX_DEPTH` to prevent infinite recursion.

            Notes:
                - The function relies on a global variable `cur_depth` to track recursion depth.
                - The function uses `getName()` to determine whether a node represents a named variable.
                - Function nodes are further traversed rather than included directly.
                - If `child_node.getName()` returns `None`, it is treated as a non-terminal and is traversed.
        """
    
        global cur_depth
        cur_depth += 1
        
        if cur_depth > cn.MAX_DEPTH:
            raise ValueError(f"Number of recursions exceeded the maximum depth, which is {cn.MAX_DEPTH}")
        
        for idx in range(ast_node.getNumChildren()):
            child_node = ast_node.getChild(idx)

            if child_node.getName() is None:
                additions = CellmlReader._identify_variables(child_node, [])
                result.extend(additions)
            else:
                if child_node.isFunction():
                    additions = CellmlReader._identify_variables(child_node, [])
                    result.extend(additions)
                else:
                    result.append(child_node.getName())
        return result



    @staticmethod
    def _get_variables(ast_node: libsbml.ASTNode) -> list[str]:
        """
            Extracts all variable names from an SBML ASTNode expression.

            This function initializes recursion depth and invokes a helper method to
            recursively traverse the AST and collect all variable names present in the expression.

            Parameters:
                ast_node (libsbml.ASTNode): The root node of the SBML Abstract Syntax Tree (AST)
                                            representing a mathematical expression.

            Returns:
                list: A list of variable names (strings) extracted from the AST.

            Notes:
                - Uses a global `cur_depth` variable to track recursion depth during traversal.
                - Delegates the recursive extraction to `_identify_variables`.
                - If the root node has a name, it is included in the results.
        """

        global cur_depth
        
        cur_depth = 0
    
        if ast_node.getName() is None:
            variables = []
        else:
            variables = [ast_node.getName()]

        result = CellmlReader._identify_variables(ast_node, variables)

        return result