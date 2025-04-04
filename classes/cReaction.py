from classes.ReactionPropertiesMixin import *


class Reaction(ReactionPropertiesMixin):

    def __init__(self, ID):

        self._ID = ID
        self._annotations = None
        self._reversible: bool = None
        self._kinetic_forward_rate_constant = None
        self._kinetic_reverse_rate_constant = None
        self._thermo_forward_rate_constant = None
        self._thermo_reverse_rate_constant = None
        self._kinetic_law = None
        self._sp_kinetic_law = None
        self._kinetic_law_type = None
        self._reactants = None
        self._products = None


    def getId(self):

        return self._ID
    
    def getReversible(self):

        return self._reversible
    
    def getKineticLaw(self):

        return self._kinetic_law
    
    def getListOfProducts(self):

        return self._products
    
    def getListOfReactants(self):

        return self._reactants