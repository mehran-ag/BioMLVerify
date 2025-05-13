from classes.cSpecies import *


class SpeciesReference(Species):

    def __init__(self, species_instance):
        # Copy all attributes from parent dynamically
        self.__dict__.update(vars(species_instance))
        self._reaction_id = None
        self._stoichiometry = None

    @property
    def stoichiometry(self):
        return self._stoichiometry
    
    @stoichiometry.setter
    def stoichiometry(self, stoichiometry):
        if isinstance(stoichiometry, (int, float)):
            self._stoichiometry = stoichiometry
        else:
            raise ValueError("Input for Stoichiometry must be a number!")
        
    @property
    def reaction_id(self):
        return self._reaction_id
    
    @reaction_id.setter
    def reaction_id(self, r_id):
        if isinstance(r_id, str):
            self._reaction_id = r_id
        else:
            raise ValueError("Input must be a string!")

    def getStoichiometry(self):

        return self._stoichiometry
    
    def getSpecies(self):

        return self._ID