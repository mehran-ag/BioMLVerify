import libsbml
import sys
import os
import matrix_constructor
import model_checker
import exceptions
import utility

from classes.cReaction import *
from classes.cModel import *
from classes.cSpecies import *
from classes.cParameter import *
from classes.cSpeciesReference import *



class BioModel(object):
    '''
    This class creates a model, either CellML or SBML, that will be processed by other classes
    '''

    def __init__(self):
        '''
        Initializes the ModelReader by reading path for a file or a container of files
        '''
        
        self._file_path = None
        self._file_name= None
        self._file_format = None
        self._biomodel = None
        self._matrix_constructor = matrix_constructor.MatrixConstructor()
        self._model_checker = model_checker.ModelChecker()


    @property
    def file_name(self):
        return self._file_name


    def read_file(self, file_path):

        try:
            if not isinstance(file_path, str):
                raise TypeError("File path must be a string.")
            
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
            
            self._file_path = file_path
            self._file_name = os.path.basename(file_path)
            self._file_format = os.path.splitext(file_path)[1][1:]
            
        except TypeError as e:
            utility.message_printer("File not read!", color="red")
            utility.printer("\nError: ", e)
            return

        except FileNotFoundError as e:
            utility.message_printer("File not read!", color="red")
            utility.printer("\nError: ", e)
            return
        
        except Exception as e:
            utility.message_printer("File not read!", color="red")
            utility.printer("\nUnexpected error: ", e)
            return

        else:
            if self._file_format == 'xml':

                utility.message_printer(f"\n\u27A4\u27A4\u27A4 The input file: {self._file_name} is a SBML model \u27A4\u27A4\u27A4", color="green", style="normal")

                self._biomodel =  self._SBML_reader()

                if self._biomodel is not None:

                    utility.message_printer(f"\n\u27A4\u27A4\u27A4 The SBML model: {self._file_name} has been succesfully converted to a BioModel\u27A4\u27A4\u27A4", color="green", style="normal")



    def _SBML_reader(self):

        """
        Reads an SBML file using libSBML.

        :return: SBML model if successful, None otherwise.
        """

        reader = libsbml.SBMLReader()
        document = reader.readSBML(self._file_path)
        if document.getNumErrors() > 0:
            utility.error_printer("\nError: ", f"The SBML file \"{self._file_name}\" contains {document.getNumErrors()} error(s).")
            utility.message_printer("\nModel not read", color="red")
            return None
        else:
            self._sbmodel = document.getModel()
            self._biomodel = Model(self._sbmodel.getId())
            self._biomodel.species = self.SBML_to_BioModel_species_tranfer(self._sbmodel)
            self._biomodel.reactions = self.SBML_to_BioModel_reaction_tranfer(self._sbmodel)
            self._biomodel.parameters = self.SBML_to_BioModel_parameter_transfer(self._sbmodel)
        
            return self._biomodel
        

    def checkModelReversibility(self, return_irreversibles = False, printing = "off"):

        try:

            if return_irreversibles:

                reversibility, irreversibles = self._model_checker.check_model_reversibility(self._biomodel, return_irreversibles = return_irreversibles)

                if printing.lower() == "on":

                    if reversibility:

                        utility.printer("\nAll reactions in the model are REVERSIBLE.\nModel: ", self._file_name)

                    else:

                        utility.printer("\nAll reactions in the model are NOT reversible.\nModel: ", self._file_name)

                return reversibility, irreversibles
            
            else:

                reversibility = self._model_checker.check_model_reversibility(self._biomodel)

                if printing.lower() == "on":

                    if reversibility:

                        utility.printer("\nAll reactions in the model are REVERSIBLE.\nModel: ", self._file_name)

                    else:

                        utility.printer("\nAll reactions in the model are NOT reversible.\nModel: ", self._file_name)

                return reversibility

        except exceptions.NoModel as e:

            utility.printer("\nAn error has been raised in function: ", "checkModelConsistency")
            utility.error_printer("Error:", e)

    def getListOfReactions(self):

        return self._biomodel._reactions
    
    def getListOfSpecies(self):

        return self._biomodel._species

        
    def getStoichiometricMatrix(self, printing = "off"):

        try:
            stoichiometric_matrix = self._matrix_constructor.stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Stoichiometric Matrix is:\n", stoichiometric_matrix)

            return stoichiometric_matrix

        except exceptions.NoModel as e:

            utility.printer("\nAn error has been raised in function: ", "getStoichiometricMatrix")
            utility.error_printer("Error:", e)

        
    def getForwardStoichiometricMatrix(self, printing = "off"):

        try:
            forward_stoichiometric_matrix = self._matrix_constructor.forward_stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Forward Stoichiometric Matrix is:\n", forward_stoichiometric_matrix)

            return forward_stoichiometric_matrix

        except exceptions.NoModel as e:

            utility.printer("\nAn error has been raised in function: ", "getForwardStoichiometricMatrix")
            utility.error_printer("Error:", e)
        

    def getReverseStoichiometricMatrix(self, printing = "off"):

        try:
            reverse_stoichiometric_matrix = self._matrix_constructor.reverse_stoichiometric_matrix_constructor(self._biomodel)

            if printing.lower() == "on":
                utility.printer("\nThe Reverse Stoichiometric Matrix is:\n", reverse_stoichiometric_matrix)

            return reverse_stoichiometric_matrix

        except exceptions.NoModel as e:

            utility.printer("\nAn error has been raised in function: ", "getReverseStoichiometricMatrix")
            utility.error_printer("Error:", e)
        
    
    def getStoichiometricColumnNamesIndices(self):

        return self._matrix_constructor.stoichiometric_matrix_column_names(self._biomodel)
    
    def getStoichiometricRowNamesIndices(self):

        return self._matrix_constructor.stoichiometric_matrix_row_names(self._biomodel)
    
    def getElementInformationInStoichiometricMatrix(self, i, j, printing = "off"):

        try:

            if self._biomodel is None:
                raise exceptions.NoModel("No BioModel has been read!!!")

            return self._matrix_constructor.stoichiometric_matrix_element_information(i, j, self._biomodel, printing = printing)
        
        except exceptions.NoModel as e:

            utility.printer("\nAn error has been raised in function: ", "getElementInformationInStoichiometricMatrix")
            utility.error_printer("Error:", e)
    
    def getThermoConversionMatrix(self, printing = "off"):

        try:
            conversion_matrix = self._matrix_constructor.kinetic_thermo_conversion_matrix_constructor(self._biomodel, printing = printing)

            return conversion_matrix
            
        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getThermoConversionMatrix")
            utility.error_printer("Error:", e)

    def getKineticRateConstantsVector(self, printing = "off"):

        try:
            kinetic_constants_vector = self._matrix_constructor.kinetic_constants_vector_constructor(self._biomodel, printing)

            return kinetic_constants_vector

        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "getKineticRateConstantsVector")
            utility.error_printer("Error:", e)

    def KineticConstantsThermoCompatibilty(self, printing = "off"):

        try:
            comatibility = self._matrix_constructor.kinetic_rates_thermo_compatibility_check(self._biomodel, printing)

            return comatibility
        
        except exceptions.NoModel as e:
            utility.printer("\nAn error has been raised in function: ", "KineticConstantsThermoCompatibilty")
            utility.error_printer("Error:", e)


    def SBML_to_BioModel_species_tranfer(self, libsbml_model):
        '''
        This function gets a SBML model, reads the required information for the species and creates a Species class for each one
        Then, it returns a list that contains the classes of species for this tool
        '''

        self._biomodel_species_list = []

        list_of_libsbml_species = libsbml_model.getListOfSpecies()



        for libsbml_species_class in list_of_libsbml_species:

            species_id = libsbml_species_class.getId()

            biomodel_species = Species(species_id)

            biomodel_species.initial_concentration = libsbml_species_class.getInitialConcentration()

            biomodel_species.compartment = libsbml_species_class.getCompartment()

            biomodel_species.charge = libsbml_species_class.getCharge()

            self._biomodel_species_list.append(biomodel_species)

        return self._biomodel_species_list
    
    def SBML_to_BioModel_parameter_transfer(self, libsbml_model):

        biomodel_parameters_list = []

        libsbml_parameters_list = libsbml_model.getListOfParameters()

        for libsbml_parameter_class in libsbml_parameters_list:

            parameter_id = libsbml_parameter_class.getId()

            parameter_value = libsbml_parameter_class.getValue()

            biomodel_parameter = Parameter(parameter_id)

            biomodel_parameter.value = parameter_value

            biomodel_parameters_list.append(biomodel_parameter)

        return biomodel_parameters_list
    

    def SBML_to_BioModel_reaction_tranfer(self, libsbml_model):

        biomodel_reactions_list = []

        libsbml_reactions = libsbml_model.getListOfReactions()

        for libsbml_reaction_class in libsbml_reactions:

            biomodel_products_list =[]
            biomodel_reactants_list = []

            libsbml_reactants = libsbml_reaction_class.getListOfReactants()

            libsbml_products = libsbml_reaction_class.getListOfProducts()

            reaction_id = libsbml_reaction_class.getId()

            index = Reaction.getCurrentIndex()

            biomodel_reaction = Reaction(reaction_id)

            biomodel_reaction.reversible = libsbml_reaction_class.getReversible()

            biomodel_reaction.kinetic_law = libsbml_reaction_class.getKineticLaw().getFormula()

            for libsbml_reactant_class in libsbml_reactants:

                id = libsbml_reactant_class.getSpecies()

                if id == "empty":

                    Reaction.ResetCounter(index)

                    biomodel_reaction.ResetIndex()

                    biomodel_reaction.boundary_condition = True

                for biomodel_species in self._biomodel_species_list:

                    if id == biomodel_species.ID:

                        biomodel_species_reference = SpeciesReference(biomodel_species)

                        biomodel_species_reference.reaction_id = reaction_id

                        biomodel_species_reference.stoichiometry = libsbml_reactant_class.getStoichiometry()

                        biomodel_reactants_list.append(biomodel_species_reference)

            for libsbml_product_class in libsbml_products:

                id = libsbml_product_class.getSpecies()

                if id == "empty":

                    Reaction.ResetCounter(index)

                    biomodel_reaction.ResetIndex()

                    biomodel_reaction.boundary_condition = True

                for biomodel_species in self._biomodel_species_list:

                    if id == biomodel_species.ID:

                        biomodel_species_reference = SpeciesReference(biomodel_species)

                        biomodel_species_reference.reaction_id = reaction_id

                        biomodel_species_reference.stoichiometry = libsbml_product_class.getStoichiometry()

                        biomodel_products_list.append(biomodel_species_reference)

            biomodel_reaction.reactants = biomodel_reactants_list

            biomodel_reaction.products = biomodel_products_list

            biomodel_reactions_list.append(biomodel_reaction)

        return biomodel_reactions_list

            



