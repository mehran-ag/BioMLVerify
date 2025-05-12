# Importing external packages
from libcellml import  Parser, Validator, Analyser, AnalyserExternalVariable, Importer, cellmlElementTypeAsString, AnalyserModel
import numpy as np
import chemparse as chp
import libchebipy as chb
import warnings

# Importing internal packages
import sys
import exceptions
import os
import utility
from pathlib import Path, PurePath

from classes.cReaction import *
from classes.cModel import *
from classes.cSpecies import *
from classes.cParameter import *
from classes.cSpeciesReference import *



class CellmlReader:

    def __init__(self):

        self._cellml_model = None
        self._model_name = None

        self._cellml_components = []
        self._cellml_variables = []

        self._variables = []
        self._coefficients = []
        self._rates = []
        self._rate_constants = []
        self._boundary_conditions = []
        self._equation_variables = []
        self._boundary_values = []
        self._enzymes = []


    def _read_file( self, file_path, cellml_strict_mode = False):

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
        
        self._model_name = os.path.basename(file_path)

        self._cellml_model = CellmlReader._read_analyse_cellml_model( file_path, cellml_strict_mode )

        self._cellml_elements_reader()

        self._variable_distinguisher()

        self._species_identifier()

        self._reaction_identifier()




    
    def _cellml_elements_reader(self):

        number_of_components = self._cellml_model.componentCount()

        if number_of_components >= 1 :

            for i in range(number_of_components):

                self._cellml_components.append( self._cellml_model.component( i ) )

        else:
            raise ValueError(f"There is no component in {self._model_name}")
        

        for component in self._cellml_components:

            number_of_variables = component.variableCount()

            for i in range(number_of_variables):

                self._cellml_variables.append(component.variable(i))

        return self._cellml_components, self._cellml_variables
    

    def _variable_distinguisher(self):

        variable_classifier = {
            'va': lambda v: self._variables.append(v),
            'co': lambda v: self._coefficients.append(v),
            'rc': lambda v: self._rate_constants.append(v),
            'ra': lambda v: self._rates.append(v),
            'bc': lambda v: self._boundary_conditions.append(v),
            'ev': lambda v: self._equation_variables.append(v),
            'bv': lambda v: self._boundary_values.append(v),
            'en': lambda v: self._enzymes.append(v)
        }

        for cellml_variable in self._cellml_variables:

            variable_id = cellml_variable.id()

            identifier = variable_id.split('_')[0]

            operation = variable_classifier.get(identifier)

            if operation:

                operation(cellml_variable)




    def _species_identifier(self):

        self._biomodel_species_list = []

        # This function determines if the characters used in a variable are all digits or alphabets
        def all_digits( word ):

            return all( char.isdigit() for char in word )

        for variable in self._variables:

            name = variable.name()

            biomodel_species = Species(name)

            name_code = variable.id().split('_')[1]

            if all_digits(name_code):

                compound, composition = CellmlReader.chebi_comp_parser(name_code)

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

            self._biomodel_species_list.append(biomodel_species)



    def _reaction_identifier(self):

        self._biomodel_reactions_list = []

        for coefficient in self._coefficients:
            
            name_code = coefficient.id().split('_')[1]

            reaction_number_parts = coefficient.id().split('_')[2]

            for reaction_number_part in reaction_number_parts.split('-'):

                reaction_number = reaction_number_part.split('.')[0]

                if not any(biomodel_reaction.ID == reaction_number for biomodel_reaction in self._biomodel_reactions_list):
                    
                    biomodel_reaction = Reaction(reaction_number)

                    self._biomodel_reactions_list.append(biomodel_reaction)



    @staticmethod
    def chebi_comp_parser( chebi_code ):

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




    @staticmethod
    def chebi_formula ( chebi_code ):
    
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




    @staticmethod
    def _read_analyse_cellml_model( file_path, cellml_strict_mode ):

        model_name = os.path.basename(file_path)
            
        parser = Parser(cellml_strict_mode)

        with open(file_path, 'r') as f:
            cellml_model = parser.parseModel(f.read())

            no_parser_warnings = parser.issueCount()

            if no_parser_warnings > 0:

                for i in range(no_parser_warnings):
                    utility.message_printer("\n" + parser.issue(i).description(), color="white", style='normal')

        validator = Validator()
        validator.validateModel(cellml_model)

        no_validator_warnings = validator.issueCount()

        if no_validator_warnings > 0:

            for i in range(no_validator_warnings):
                utility.message_printer("\n" + validator.issue(i).description(), color="magenta", style="normal")

            raise ValueError(f"Model {model_name} has validation issues and cannot be imported!")


        importer = Importer(cellml_strict_mode)

        base_dir = os.path.dirname(file_path)
        importer.resolveImports(cellml_model, base_dir)

        no_importer_warnings = importer.issueCount()

        if no_importer_warnings > 0:

            for i in range(no_importer_warnings):
                utility.message_printer("\n" + importer.issue(i).description(), color="magenta", style="normal")

            raise ValueError(f"Model {model_name} has import issues and cannot be imported!")

        else:
            if cellml_model.hasUnresolvedImports():
                utility.message_printer("There are Unresolved Import Issues", color="magenta", style="normal")

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
                utility.message_printer("\n" + analyser.issue(i).description(), color="magenta", style="normal")

            #raise ValueError(f"Model {model_name} has analyse issues and cannot be imported!")
        

        return cellml_model



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