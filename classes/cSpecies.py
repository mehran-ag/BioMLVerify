from classes.SpeciesPropertiesMixin import *


class Species(SpeciesPropertiesMixin):

    _counter: int = 0

    def __init__(self, ID):

        if ID != "empty":
            self._index = Species._counter
            Species._counter += 1
        else:
            self._index = None

        self._ID = ID
        self._name = None
        self._initial_concentration = None
        self._compartment = None
        self._annotations = None
        self._charge = None
        self._thermodynamic_rate_constant = None
        self._compound = None  #This is the scientific name of the species like CO2, H2O, CH4
        self._composition = None
        self._chebi_code = None

    @classmethod
    def getCurrentIndex(cls):
        return cls._counter

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
    

    @staticmethod
    def reset_counter():
        Species._counter = 0