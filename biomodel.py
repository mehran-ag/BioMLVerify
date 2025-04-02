import libsbml
import sys
import os
from colorama import Fore, Back, Style, init
import matrix_constructor as mc
import exceptions

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
        
        self.file_path = None
        self.file_name= None
        self.file_format = None
        self.model = None
        self.matrix_constructor = mc.MatrixConstructor()


    def read_file(self, file_path):
                
        init( autoreset=True )

        try:
            if not isinstance(file_path, str):
                raise TypeError(Fore.RED + "File path must be a string.")
            
            if not os.path.exists(file_path):
                raise FileNotFoundError(Fore.RED + f"File not found: {file_path}")
            
            self.file_path = file_path
            self.file_name = os.path.basename(file_path)
            self.file_format = os.path.splitext(file_path)[1][1:]
            
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
            if self.file_format == 'xml':

                self.model =  self._SBML_reader()

                print(Fore.GREEN + "\n\u27A4\u27A4\u27A4 The input file is a SBML model \u27A4\u27A4\u27A4")



    def _SBML_reader(self):

        """
        Reads an SBML file using libSBML.

        :return: SBML model if successful, None otherwise.
        """

        reader = libsbml.SBMLReader()
        document = reader.readSBML(self.file_path)
        if document.getNumErrors() > 0:
            print(f"\nError: The SBML file contains {document.getNumErrors()} error(s).")
            print("\nModel not read")
            return
        else:
            self.sbmodel = document.getModel()
            biomodel = Model(self.sbmodel.getId())
            biomodel.species = self.SBML_to_BioModel_species_tranfer(self.sbmodel)
            biomodel.reactions = self.SBML_to_BioModel_reaction_tranfer(self.sbmodel)
            biomodel.parameters = self.SBML_to_BioModel_parameter_transfer(self.sbmodel)
            
            return biomodel
            #return self.sbmodel
        
    def getStoichiometricMatrix(self):

        try:
            stoichiometic_matrix = self.matrix_constructor.SBML_stoichiomrtic_matrix_constructor(self.model)

            return stoichiometic_matrix

        except exceptions.NoModel as e:

            print(Fore.RED + f"\nError: {e}")

            return "Error encountered"
        
    
    def getStoichiometricColumnNamesIndices(self):

        return self.matrix_constructor.stoichiometric_matrix_column_names()
    
    def getStoichiometricRowNamesIndices(self):

        return self.matrix_constructor.stoichiometric_matrix_row_names()
    
    def getElementInformationInStoichiometricMatrix(self, i,j):

        return self.matrix_constructor.elementinformation(i,j)
    
    def getThermoConversionMatrix(self):

        try:
            self.matrix_constructor.conversion_matrix_constructor(self.model)
            
        except exceptions.NoModel as e:
            print(Fore.BLUE + "\nAn error has been raised in \"getThermoConversionMatrix\" function")
            print(Fore.RED + f"{e}")
    


    def SBML_to_BioModel_species_tranfer(self, model):

        self.biomodel_species_list = []

        list_of_species = model.getListOfSpecies()

        for species_class in list_of_species:

            species_id = species_class.getId()

            biomodel_species = Species(species_id)

            biomodel_species.initial_concentration = species_class.getInitialConcentration()

            biomodel_species.compartment = species_class.getCompartment()

            biomodel_species.charge = species_class.getCharge()

            self.biomodel_species_list.append(biomodel_species)

        return self.biomodel_species_list
    
    def SBML_to_BioModel_parameter_transfer(self, model):

        biomodel_parameters_list = []

        parameters_list = model.getListOfParameters()

        for parameter_class in parameters_list:

            parameter_id = parameter_class.getId()

            biomodel_parameter = Parameter(parameter_id)

            biomodel_parameters_list.append(biomodel_parameter)

        return biomodel_parameters_list
    

    def SBML_to_BioModel_reaction_tranfer(self, model):

        biomodel_reactions_list = []

        reactions = model.getListOfReactions()

        for reaction_class in reactions:

            biomodel_products_list =[]
            biomodel_reactants_list = []

            reactants = reaction_class.getListOfReactants()

            products = reaction_class.getListOfProducts()

            reaction_id = reaction_class.getId()

            biomodel_reaction = Reaction(reaction_id)

            biomodel_reaction.reversible = reaction_class.getReversible()

            biomodel_reaction.kinetic_law = reaction_class.getKineticLaw().getFormula()

            for reactant_class in reactants:

                id = reactant_class.getSpecies()

                for biomodel_species in self.biomodel_species_list:

                    if id == biomodel_species.ID:

                        species_reference = SpeciesReference(biomodel_species)

                        species_reference.reaction_id = reaction_id

                        species_reference.stoichiometry = reactant_class.getStoichiometry()

                        biomodel_reactants_list.append(species_reference)

            for product_class in products:

                id = product_class.getSpecies()

                for biomodel_species in self.biomodel_species_list:

                    if id == biomodel_species.ID:

                        species_reference.reaction_id = reaction_id

                        species_reference.stoichiometry = reactant_class.getStoichiometry()

                        biomodel_products_list.append(species_reference)

            biomodel_reaction.reactants = biomodel_reactants_list

            biomodel_reaction.products = biomodel_products_list

            biomodel_reactions_list.append(biomodel_reaction)

        return biomodel_reactions_list

            



