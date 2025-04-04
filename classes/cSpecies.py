from classes.SpeciesPropertiesMixin import *


class Species(SpeciesPropertiesMixin):

    def __init__(self, ID):
        self._ID = ID
        self._name = None
        self._initial_concentration = None
        self._compartment = None
        self._annotations = None
        self._charge = None

    def getId(self):

        return self._ID
    
    def getName(self):

        return self._ID
    
    def getInitialConcentration(self):

        return self._initial_concentration
    
    def getCharge(self):

        return self._charge
    
    def getCompartment(self):

        return self._compartment