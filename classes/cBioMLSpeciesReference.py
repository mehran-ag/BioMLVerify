from classes.cBioMLSpecies import *


class BioMLSpeciesReference(BioMLSpecies):

    def __init__(self, species_instance):
        # Copy all attributes from parent dynamically
        self.__dict__.update(vars(species_instance))
        self._reaction_id: str = None
        self._stoichiometry: float = None

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
        

        

    def get_stoichiometry(self):

        return self._stoichiometry
    
    def get_species(self):

        return self._ID