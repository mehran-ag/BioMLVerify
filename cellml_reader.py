# Importing external packages
from libcellml import  Parser, Validator, Analyser, AnalyserExternalVariable, Importer, cellmlElementTypeAsString, AnalyserModel
import numpy as np
import chemparse as chp
import warnings

# Importing internal packages
import sys
import exceptions
import os
import utility
from pathlib import Path, PurePath


class CellmlReader:


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
        
        model_name = os.path.basename(file_path)

        cellml_model = _read_analyse_cellml_model( file_path, cellml_strict_mode )

        number_of_components = cellml_model.componentCount()

        if number_of_components >= 1 :

            components = []

            for component in range(0, number_of_components):

                components.append( cellml_model.component( component ) )

        else:
            print( 'There is no component in the CellML file' )
            exit(-4)










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

        raise ValueError(f"Nodel {model_name} has import issues and cannot be imported!")

    else:
        if cellml_model.hasUnresolvedImports():
            utility.message_printer("There are Unresolved Import Issues", color="magenta", style="normal")

            raise ValueError(f"Nodel {model_name} has import issues and cannot be imported!")


    flat_cellml_model = importer.flattenModel(cellml_model)

    if not flat_cellml_model:
        raise ValueError(f"Model {model_name} cannot be imported: flattening issues!")
    
    external_variables_dic = _ext_var_dic(flat_cellml_model)

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

        raise ValueError(f"Nodel {model_name} has analyse issues and cannot be imported!")
    

    return cellml_model




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