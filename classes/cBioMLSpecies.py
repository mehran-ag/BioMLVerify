from classes.BioMLSpeciesPropertiesMixin import *


class BioMLSpecies(BioMLSpeciesPropertiesMixin):

    _counter: int = 0

    def __init__(self, ID):

        if ID != "empty":
            self._index = BioMLSpecies._counter
            BioMLSpecies._counter += 1
        else:
            self._index: int = None

        self._ID: str = ID
        self._name: str = None
        self._initial_concentration: float = None
        self._compartment: str = None
        self._annotations: dict[str, list[str]] = {}
        self._charge: int = 0
        self._thermodynamic_rate_constant: str = None
        self._compound: str = None  #This is the scientific name of the species like CO2, H2O, CH4
        self._composition: dict = None
        self._chebi_code: str = None

    @classmethod
    def get_current_index(cls):
        return cls._counter

    def get_id(self):

        return self._ID
    
    def get_name(self):

        return self._ID
    
    def get_initial_concentration(self):

        return self._initial_concentration
    
    def get_charge(self):

        return self._charge
    
    def get_compartment(self):

        return self._compartment
    
    @classmethod
    def reset_counter(cls, new_counter_value = 0):
        cls._counter = new_counter_value  # Reset the class-level counter