# Importing external packages
from libcellml import  Parser, Validator, Analyser, AnalyserExternalVariable, Importer, cellmlElementTypeAsString, AnalyserModel
import numpy as np
import chemparse as chp
import libchebipy as chb
import libsbml
import sympy as sp

from lxml import etree

# Importing internal packages
import sys
import exceptions
import os
import utility
from pathlib import Path, PurePath
import constants as cn

from xml.dom.minidom import parseString

from classes.cReaction import *
from classes.cModel import *
from classes.cSpecies import *
from classes.cParameter import *
from classes.cSpeciesReference import *



class CellmlReader:



    # ********************************
    # *           Function           *
    # ********************************
    def read_file( self, file_path, cellml_strict_mode = False):

        base_dir = os.path.dirname(file_path)

        if not isinstance(base_dir, str):
            if isinstance(base_dir, Path):
                base_dir = str(base_dir)
            elif isinstance(base_dir,PurePath):
                base_dir = str(base_dir)
            else:
                raise TypeError('The file path should be a string or a Path object.')

        if not os.path.isfile(file_path):
            raise FileNotFoundError('Model source file `{}` does not exist.'.format(file_path))
        
        cellml_model_name = os.path.basename(file_path)

        cellml_model = CellmlReader._read_analyse_cellml_model( file_path, cellml_strict_mode )

        cellml_contents = self._extract_cellml_content(cellml_model, cellml_model_name) #returns a dictionary containing "cellml_eqs", "cellml_species_instances", "cellml_vars_instances", and "cellml_ast_nodes"

        cellml_vars_instances = cellml_contents["cellml_vars_instances"]

        variable_type_buckets = self._variable_distinguisher(cellml_vars_instances)
        
        '''
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
        '''

        keys_to_check = ['va', 'co', 'rc', 'ra']

        if all(variable_type_buckets.get(k) for k in keys_to_check):


            biomodel_species_list = self._species_identifier(variable_type_buckets['va'])

            biomodel_reactions_list = self._reaction_identifier(variable_type_buckets['co'], biomodel_species_list)

            # biomodel_reactions_list, biomodel_species_list = self._boundary_condition_identifier(variable_type_buckets['bc'], biomodel_reactions_list, biomodel_species_list)

            self._kinetic_rate_constant_finder(variable_type_buckets['rc'], biomodel_reactions_list)

            biomodel = Model(cellml_model.id())

            biomodel.species = biomodel_species_list

            biomodel.reactions = biomodel_reactions_list

            mass_action = self._find_cellml_mass_actions(**cellml_contents)

            if mass_action:
            
                biomodel.is_mass_action = True

            else:

                biomodel.is_mass_action = False

            return biomodel
            
        else:

            mass_action = self._find_cellml_mass_actions(**cellml_contents)

            if mass_action:

                print("It is mass action")

            else:

                print(f"\nCellML model \"{Path(cellml_model_name).stem.upper()}\" has (an) equation(s) not governed by Mass Action Kinetics\n")

                return

        



    # ********************************
    # *           Function           *
    # ********************************
    def read_file_old( self, file_path, cellml_strict_mode = False):

        base_dir = os.path.dirname(file_path)

        if not isinstance(base_dir, str):
            if isinstance(base_dir, Path):
                base_dir = str(base_dir)
            elif isinstance(base_dir,PurePath):
                base_dir = str(base_dir)
            else:
                raise TypeError('The file path should be a string or a Path object.')

        if not os.path.isfile(file_path):
            raise FileNotFoundError('Model source file `{}` does not exist.'.format(file_path))
        
        cellml_model_name = os.path.basename(file_path)

        cellml_model = CellmlReader._read_analyse_cellml_model( file_path, cellml_strict_mode )

        cellml_components, cellml_variables = self._cellml_elements_reader(cellml_model, cellml_model_name)

        variable_type_buckets = self._variable_distinguisher(cellml_variables)
        
        '''
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
        '''

        biomodel_species_list = self._species_identifier(variable_type_buckets['va'])

        biomodel_reactions_list = self._reaction_identifier(variable_type_buckets['co'], biomodel_species_list)

        # biomodel_reactions_list, biomodel_species_list = self._boundary_condition_identifier(variable_type_buckets['bc'], biomodel_reactions_list, biomodel_species_list)

        self._kinetic_rate_constant_finder(variable_type_buckets['rc'], biomodel_reactions_list)

        biomodel = Model(cellml_model.id())

        biomodel.species = biomodel_species_list

        biomodel.reactions = biomodel_reactions_list

        return biomodel




    # ********************************
    # *           Function           *
    # ********************************
    def _cellml_elements_reader(self, cellml_model, cellml_model_name):
        '''
        
        '''

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
    def _variable_distinguisher(self, cellml_variables):

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
    def _species_identifier(self, variables):

        biomodel_species_list = []

        for variable in variables:

            name = variable.name()

            matched_species = next( ( species_instance for species_instance in biomodel_species_list if species_instance.ID == name ), None )

            if matched_species is not None:
                continue

            biomodel_species = Species(name)

            name_code = variable.id().split('_')[1]

            if all( char.isdigit() for char in name_code ):

                compound, composition = CellmlReader._chebi_comp_parser(name_code)

                biomodel_species.chebi_code = name_code

            else:

                # If there is no chebi code for the variable, then we have to split the code more
                # A '-' after the name of the compound shows its composition, and each species is separated by a '.'
                if ( len(name_code.split('-')) == 1 ):
                    compound = name_code
                    composition={ name_code: 1 }
                
                else:

                    compound = name_code.split('-')[0]
                    formula = name_code.split('-')[1].split('.')

                    composition={}
                    for element in formula:
                        composition[element] = 1    #IT SHOULD BE NOTED THAT CURRENTLY I CAN INCLUDE ONLY ONE INSTANCE OF EACH ELEMENT IN THE COMPOSITION. IT IS THE DEFAULT VALUE SET HERE
                                                    #ATTENTION ATTENTION

            biomodel_species.compound = compound
            biomodel_species.composition = composition

            biomodel_species_list.append(biomodel_species)

        return biomodel_species_list





    # ********************************
    # *           Function           *
    # ********************************
    def _reaction_identifier(self, coefficients, biomodel_species_list ):

        biomodel_reactions_list = []

        for coefficient in coefficients:
            
            name_code = coefficient.id().split('_')[1]

            if all( char.isdigit() for char in name_code ):

                compound, _ = CellmlReader._chebi_comp_parser(name_code)

            else:

                compound = name_code

            matched_biomodel_species =  next( ( species_instance for species_instance in biomodel_species_list if species_instance.compound == compound ), None )

            if matched_biomodel_species is None:
                raise ValueError(f"There is no match for species {compound} in the list of species")
            

            reaction_number_parts = coefficient.id().split('_')[2]

            for i, reaction_number_part in enumerate(reaction_number_parts.split('-')):

                reaction_number = reaction_number_part.split('.')[0]

                matched_biomodel_reaction = next( ( reaction_instance for reaction_instance in biomodel_reactions_list if reaction_instance.ID == reaction_number ), None )

                if matched_biomodel_reaction is None:
                    
                    biomodel_reaction = Reaction(reaction_number)

                    biomodel_species_ref = SpeciesReference(matched_biomodel_species)

                    biomodel_species_ref.reaction_id = reaction_number

                    stoichiometry = float(coefficient.initialValue())

                    if stoichiometry > 0:

                        biomodel_species_ref.stoichiometry = abs(stoichiometry)

                        biomodel_reaction.products.append(biomodel_species_ref)

                    elif stoichiometry < 0:

                        biomodel_species_ref.stoichiometry = abs(stoichiometry)

                        biomodel_reaction.reactants.append(biomodel_species_ref)

                    else:

                        raise ValueError(f"The value of stoichiometric coefficient is not acceptable: {coefficient.initialValue()}")

                    biomodel_reactions_list.append(biomodel_reaction)

                else:

                    if i > 0:
                        continue

                    biomodel_species_ref = SpeciesReference(matched_biomodel_species)

                    biomodel_species_ref.reaction_id = reaction_number

                    stoichiometry = float(coefficient.initialValue())

                    if stoichiometry > 0:

                        biomodel_species_ref.stoichiometry = abs(stoichiometry)

                        matched_biomodel_reaction.products.append(biomodel_species_ref)

                    elif stoichiometry < 0:

                        biomodel_species_ref.stoichiometry = abs(stoichiometry)

                        matched_biomodel_reaction.reactants.append(biomodel_species_ref)

                    else:

                        raise ValueError(f"The value of stoichiometric coefficient is not acceptable: {coefficient.initialValue()}")
                    
        return biomodel_reactions_list







    # ********************************
    # *           Function           *
    # ********************************
    def _boundary_condition_identifier(self, boundary_conditions, biomodel_reactions_list, biomodel_species_list):

        if boundary_conditions:

            for bc in boundary_conditions:

                reaction_number = bc.id()

                if not any(biomodel_reaction.ID == reaction_number for biomodel_reaction in biomodel_reactions_list):
                        
                    biomodel_reaction = Reaction(reaction_number)

                    biomodel_reaction.boundary_condition = True

                    name_code =  bc.id().split('_')[1].split('.')[0]

                    try:

                        flow_direction = bc.id().split('_')[2].split('.')[0]

                    except IndexError:

                        flow_direction = ''

                    if all( char.isdigit() for char in name_code ):

                        compound, _ = CellmlReader._chebi_comp_parser(name_code)

                    else:

                        compound = name_code


                    matched_species = next( ( species_instance for species_instance in biomodel_species_list if species_instance.compound == compound ), None )

                    if matched_species is None:
                        raise ValueError(f"The Species for the boundary condition {reaction_number} can't be found!")
                    
                    
                    biomodel_species_ref = SpeciesReference(matched_species)

                    biomodel_species_ref.reaction_id = reaction_number

                    biomodel_species_ref.stoichiometry = 1

                    if flow_direction == 'i':

                        biomodel_reaction.products.append(biomodel_species_ref)

                    elif flow_direction == 'o':

                        biomodel_reaction.reactants.append(biomodel_species_ref)

                    else:

                        biomodel_reaction.reactants.append(biomodel_species_ref)
                    

                    # Now, we will create the virtual internal and external species for the compound
                    compound_bc = compound + '_e'

                    matched_bc_species =  next( ( species_instance for species_instance in self._biomodel_species_list if species_instance.compound == compound_bc ), None )

                    if matched_bc_species is None:

                        biomodel_species = Species(compound_bc)

                        biomodel_species.compound = compound_bc

                        biomodel_species_list.append(biomodel_species)

                    else:

                        biomodel_species = matched_bc_species

                    biomodel_species_ref = SpeciesReference(biomodel_species)

                    biomodel_species_ref.reaction_id = reaction_number

                    biomodel_species_ref.stoichiometry = 1

                    if flow_direction == 'i':

                        biomodel_reaction.reactants.append(biomodel_species_ref)

                    elif flow_direction == 'o':

                        biomodel_reaction.products.append(biomodel_species_ref)

                    else:

                        biomodel_reaction.products.append(biomodel_species_ref)

                    biomodel_reactions_list.append(biomodel_reaction)

        return biomodel_reactions_list, biomodel_species_list






    # ********************************
    # *           Function           *
    # ********************************
    def _kinetic_rate_constant_finder(self, rate_constants, biomodel_reactions_list):

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

                matched_reaction =  next( ( reaction_instance for reaction_instance in biomodel_reactions_list if reaction_instance.ID == reaction_name ), None )

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
    def _extract_cellml_content(self, cellml_model, cellml_model_name):

        cellml_vars_instances = []

        cellml_species_instances = []

        cellml_eqs = []

        cellml_ast_nodes = []

        for i in range(cellml_model.componentCount()):

            component = cellml_model.component(i)

            component_name = component.name()

            component_id = component.id()

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

        return {
            "cellml_eqs": cellml_eqs,
            "cellml_species_instances": cellml_species_instances,
            "cellml_vars_instances": cellml_vars_instances,
            "cellml_ast_nodes": cellml_ast_nodes
        }




    def _find_cellml_mass_actions(self, **contents):

        cellml_species_instances = contents["cellml_species_instances"]

        cellml_vars_instances = contents["cellml_vars_instances"]

        cellml_ast_nodes = contents["cellml_ast_nodes"]

        species_in_cellml_eq = []

        mass_action = True

        cellml_vars = [var_instance.name() for var_instance in cellml_vars_instances]

        cellml_species = [species_instance.name() for species_instance in cellml_species_instances]

        # Prepare a space-separated and comma-separated version
        symbol_dict = {cellml_var: sp.symbols(cellml_var) for cellml_var in cellml_vars}

        for cellml_ast_node in cellml_ast_nodes:

            cellml_eq_vars = CellmlReader._get_variables(cellml_ast_node)

            for cellml_eq_var in cellml_eq_vars:

                if cellml_eq_var in cellml_species:

                    species_in_cellml_eq.append(cellml_eq_var)

            cellml_eq = libsbml.formulaToL3String(cellml_ast_node)

            if '==' in cellml_eq:

                cellml_eq_rhs = cellml_eq.split('==')[1].strip()


            symp_cellml_eq = sp.sympify(cellml_eq_rhs.replace("^", "**"), locals = symbol_dict)

            simp_cellml_eq = str(sp.simplify(symp_cellml_eq))

            eq_mass_action = self._check_mass_action( cellml_eq_rhs, simp_cellml_eq, cellml_vars, species_in_cellml_eq )

            if not eq_mass_action:

                mass_action = eq_mass_action

        return mass_action


    


    def _check_mass_action(self, cellml_eq, simp_cellml_eq, cellml_vars, species_in_cellml_eq):

        flag = False

        if self._single_product( cellml_eq, simp_cellml_eq ) or \
            self._diff_of_products( cellml_eq, simp_cellml_eq ):

                if len( species_in_cellml_eq ) != 0:

                    flag = True

                    if len( species_in_cellml_eq ) == 1:

                        if cellml_eq.count(species_in_cellml_eq[0]) != 1 or simp_cellml_eq.count(species_in_cellml_eq[0]) != 1:
                            flag = False
        
        try:
            fracs = self._frac_parts(simp_cellml_eq, cellml_vars)
            if len(species_in_cellml_eq) > 0:
                for i in range(len(species_in_cellml_eq)):
                    if species_in_cellml_eq[i] in fracs["denominator"]:
                        flag = False

        except:
            pass

        return flag





    def _single_product(self, kinetic_law, simple_kinetic_law):

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
    def _diff_of_products(self, kinetic_law, simple_kinetic_law):

        flag = False

        if self._single_product(kinetic_law, simple_kinetic_law) == False:
            terms = kinetic_law.split("-")
            if len(terms) == 2:
                flag = True

        return flag
    

    # ********************************
    # *           Function           *
    # ********************************
    def _frac_parts( self, simp_cellml_eq, cellml_vars ):

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
    def _has_multiple_equations(self, mathml_str):
        """
        Checks whether there are more than one <apply><eq/>...</apply> elements in the MathML.
        """
        root = etree.fromstring(mathml_str.encode())
        NSMAP = {'m': 'http://www.w3.org/1998/Math/MathML'}
        eq_applies = root.xpath(".//m:apply[m:eq]", namespaces=NSMAP)
        return len(eq_applies) > 1
    





    # ********************************
    # *           Function           *
    # ********************************
    def _split_equations(self, mathml_str):
        """
        Splits a MathML script containing multiple equations into individual <math> blocks.
        Returns a list of MathML strings.
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
    def _remove_units_from_mathml(seflf, mathml_str):
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
    @staticmethod
    def _chebi_comp_parser( chebi_code ):

        """
        This function receives a string which includes ChEBI code for the compound and uses EBI API to search for the compounds chemical composition and fetches it
        Then using the chemparse package decomposes the chemical formula to its elements and returns it as a dictionary
        It should be mentioned that some ChEBI compounds don't have chemical formula registered for them and some have two different chemical compositions assigned to them
        """

        chebi_entity = chb.ChebiEntity(chebi_code)
        name = chebi_entity.get_name()

        parsed_compound={}


        formulae = chebi_entity.get_formulae() # There are two methods to fetch chemical composition of a compound: get_formulae is used when more than one composition is registerd for the compound
                                                # It returns a list containing 'Formula' classes and Getting list's length will show us that if there is one or more chemical compositions registered
        if len(formulae) == 0: # Sometimes no chemical formula is registered so the length of the list will be zero

            parsed_compound = None
            formula = None

        elif len(formulae) == 1:

            formula = chebi_entity.get_formula()
            parsed_compound = chp.parse_formula(formula)

        else:

            formula = formulae[1].get_formula()
            parsed_compound = chp.parse_formula(formula)

        return formula, parsed_compound






    # ********************************
    # *           Function           *
    # ********************************
    @staticmethod
    def _chebi_formula ( chebi_code ):
    
        """
        This function receives a string which includes ChEBI code for the compound and uses EBI API to search for the compound's chemical composition and fetches it
        It should be mentioned that some ChEBI compounds don't have chemical formula registered for them and some have two different chemical compositions assigned to them
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





    # ********************************
    # *           Function           *
    # ********************************
    @staticmethod
    def _read_analyse_cellml_model( file_path, cellml_strict_mode ):

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

            #raise ValueError(f"Model {model_name} has analyse issues and cannot be imported!")
        

        return cellml_model








    # ********************************
    # *           Function           *
    # ********************************
    @staticmethod
    def _ext_var_dic(flatModel,external_variables_info={}):
        """
        Create a dictionary of external variables in the flattened model.
        
        Parameters
        ----------
        flatModel: Model
            The flattened CellML model.
        external_variables_info: dict
            The external variables to be specified, in the format of {id:{'component': , 'name': }}

        Raises
        ------
        ValueError
            If an external variable is not found in the flattened model.

        Returns
        -------
        dict
            The dictionary of external variables in the flattened model, in the format of {external_variable:[]}
        
        Notes
        -----
            No dependency is specified for the external variables.
        """

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
    def _variable_identifier(ast_node, result):
    
        global cur_depth
        cur_depth += 1
        
        if cur_depth > cn.MAX_DEPTH:
            raise ValueError(f"Number of recursions exceeded the maximum depth, which is {cn.MAX_DEPTH}")
        
        for idx in range(ast_node.getNumChildren()):
            child_node = ast_node.getChild(idx)

            if child_node.getName() is None:
                additions = CellmlReader._variable_identifier(child_node, [])
                result.extend(additions)
            else:
                if child_node.isFunction():
                    additions = CellmlReader._variable_identifier(child_node, [])
                    result.extend(additions)
                else:
                    result.append(child_node.getName())
        return result



    @staticmethod
    def _get_variables(ast_node):


        global cur_depth
        
        cur_depth = 0
    
        if ast_node.getName() is None:
            variables = []
        else:
            variables = [ast_node.getName()]

        result = CellmlReader._variable_identifier(ast_node, variables)

        return result