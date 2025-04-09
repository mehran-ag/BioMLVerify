import libsbml
import sys
import os
from colorama import Fore, Back, Style, init
import matrix_constructor as mc
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
        self._bio_model = None
        self._matrix_constructor = mc.MatrixConstructor()


    def read_file(self, file_path):
                
        init( autoreset=True )

        try:
            if not isinstance(file_path, str):
                raise TypeError(Fore.RED + "File path must be a string.")
            
            if not os.path.exists(file_path):
                raise FileNotFoundError(Fore.RED + f"File not found: {file_path}")
            
            self._file_path = file_path
            self._file_name = os.path.basename(file_path)
            self._file_format = os.path.splitext(file_path)[1][1:]
            
        except TypeError as e:
            print(Fore.RED + "File not read!")
            print(Fore.RED + f"Error: {e}")
            return

        except FileNotFoundError as e:
            print(Fore.RED + "File not read!")
            print(Fore.RED + f"Error: {e}")
            return
        
        except Exception as e:
            print(Fore.RED + "File not read!")
            print(Fore.RED + f"Unexpected error: {e}")
            return

        else:
            if self._file_format == 'xml':

                print(Fore.GREEN + "\n\u27A4\u27A4\u27A4 The input file is a SBML model \u27A4\u27A4\u27A4")

                self._bio_model =  self._SBML_reader()

                print(Fore.GREEN + "\n\u27A4\u27A4\u27A4 The SBML model has been succesfully converted to a BioModel\u27A4\u27A4\u27A4")



    def _SBML_reader(self):

        """
        Reads an SBML file using libSBML.

        :return: SBML model if successful, None otherwise.
        """

        reader = libsbml.SBMLReader()
        document = reader.readSBML(self._file_path)
        if document.getNumErrors() > 0:
            print(f"\nError: The SBML file contains {document.getNumErrors()} error(s).")
            print("\nModel not read")
            return
        else:
            self._sbmodel = document.getModel()
            self._biomodel = Model(self._sbmodel.getId())
            self._biomodel.species = self.SBML_to_BioModel_species_tranfer(self._sbmodel)
            self._biomodel.reactions = self.SBML_to_BioModel_reaction_tranfer(self._sbmodel)
            self._biomodel.parameters = self.SBML_to_BioModel_parameter_transfer(self._sbmodel)
        
            return self._biomodel
            #return self.sbmodel

    def getListOfReactions(self):

        return self._biomodel._reactions
    
    def getListOfSpecies(self):

        return self._biomodel._species

        
    def getStoichiometricMatrix(self, printing = "off"):

        try:
            stoichiometric_matrix = self._matrix_constructor.stoichiometric_matrix_constructor(self._bio_model)

            if printing.lower() == "on":
                utility.printer("\nThe Stoichiometric Matrix is:\n", stoichiometric_matrix)

            return stoichiometric_matrix

        except exceptions.NoModel as e:

            print(Fore.RED + f"\nError: {e}")

            return "Error encountered"
        
    def getForwardStoichiometricMatrix(self, printing = "off"):

        try:
            forward_stoichiometric_matrix = self._matrix_constructor.forward_stoichiometric_matrix_constructor(self._bio_model)

            if printing.lower() == "on":
                utility.printer("\nThe Forward Stoichiometric Matrix is:\n", forward_stoichiometric_matrix)

            return forward_stoichiometric_matrix

        except exceptions.NoModel as e:

            print(Fore.RED + f"\nError: {e}")

            return "Error encountered"
        

    def getReverseStoichiometricMatrix(self, printing = "off"):

        try:
            reverse_stoichiometric_matrix = self._matrix_constructor.reverse_stoichiometric_matrix_constructor(self._bio_model)

            if printing.lower() == "on":
                utility.printer("\nThe Reverse Stoichiometric Matrix is:\n", reverse_stoichiometric_matrix)

            return reverse_stoichiometric_matrix

        except exceptions.NoModel as e:

            print(Fore.RED + f"\nError: {e}")

            return "Error encountered"
        
    
    def getStoichiometricColumnNamesIndices(self):

        return self._matrix_constructor.stoichiometric_matrix_column_names(self._biomodel)
    
    def getStoichiometricRowNamesIndices(self):

        return self._matrix_constructor.stoichiometric_matrix_row_names(self._biomodel)
    
    def getElementInformationInStoichiometricMatrix(self, i, j):

        return self._matrix_constructor.stoichiometric_matrix_element_information(i, j, self._biomodel)
    
    def getThermoConversionMatrix(self, printing = "off"):

        try:
            conversion_matrix = self._matrix_constructor.kinetic_thermo_conversion_matrix_constructor(self._bio_model, printing = printing)

            return conversion_matrix
            
        except exceptions.NoModel as e:
            print(Fore.BLUE + "\nAn error has been raised in \"getThermoConversionMatrix\" function")
            print(Fore.RED + f"{e}")

    def getKineticRateConstantsVector(self, printing = "off"):

        try:
            kinetic_constants_vector = self._matrix_constructor.kinetic_constants_vector_constructor(self._bio_model, printing)

            return kinetic_constants_vector

        except exceptions.NoModel as e:
            print(Fore.BLUE + "\nAn error has been raised in \"ggetKineticRateConstantsVector\" function")
            print(Fore.RED + f"{e}")

    def KineticConstantsThermoCompatibilty(self, printing = "off"):

        try:
            comatibility = self._matrix_constructor.kinetic_rates_thermo_compatibility_check(self._bio_model, printing)

            return comatibility
        
        except exceptions.NoModel as e:
            print(Fore.BLUE + "\nAn error has been raised in \"KineticConstantsThermoCompatibilty\" function")
            print(Fore.RED + f"{e}")


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

            



