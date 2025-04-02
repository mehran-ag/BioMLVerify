from classes.cSpecies import *


class SpeciesReference(Species):

    def __init__(self, species_instance):
        # Copy all attributes from parent dynamically
        self.__dict__.update(vars(species_instance))
        self.reaction_id = None
        self.stoichiometry = None

    def getStoichiometry(self):

        return self.stoichiometry
    
    def getSpecies(self):

        return self.ID