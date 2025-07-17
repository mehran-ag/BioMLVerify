from classes.BioMLReactionPropertiesMixin import *
from typing import Union


class BioMLReaction(BioMLReactionPropertiesMixin):

    _counter: int = 0 #Class variable to track the index

    def __init__(self, ID):

        self._index = BioMLReaction._counter
        BioMLReaction._counter += 1
        self._ID: str = ID
        self._annotations: list = None
        self._reversible: bool = None
        self._kinetic_forward_rate_constant: str = None
        self._kinetic_forward_rate_constant_value: Union[int, float] = None
        self._kinetic_reverse_rate_constant: str = None
        self._kinetic_reverse_rate_constant_value: Union[int, float] = None
        self._thermo_forward_rate_constant: str = None
        self._thermo_forward_rate_constant_value: Union[int, float] = None
        self._thermo_reverse_rate_constant: str = None
        self._thermo_reverse_rate_constant_value: Union[int, float] = None
        self._kappa: Union[int, float] = None
        self._kinetic_law: str = None
        self._sp_kinetic_law = None #sympy expression
        self._expanded_kinetic_law: str = None #Not set yet
        self._kinetic_law_type: str = None
        self._reactants: list = []
        self._products: list = []
        self._boundary_condition: bool = False
        self._local_parameters: list = None
        self._klaw_variables: list = []
        self._mass_action: bool = None


    @classmethod
    def get_current_index(cls):
        return cls._counter
    
    @classmethod
    def reset_counter(cls, new_counter_value = 0):
        cls._counter = new_counter_value  # Reset the class-level counter




    def reset_index(self):
        self._index = None


    def assign_index(self):
        if self._index != None:
            print(f"This reaction has analready has an index: {self._index}")
        else:
            self._index = BioMLReaction._counter
            BioMLReaction._counter += 1
            print(f"Index \"{self._index}\" has now been assigned to this reaction")




    def get_id(self):

        return self._ID
    
    def get_kinetic_forward_rate_constant_value(self):
        
        return self._kinetic_forward_rate_constant_value
    
    def get_kinetic_forward_rate_constant(self):
        
        return self._kinetic_forward_rate_constant

    def get_reversible(self):

        return self._reversible
    
    def get_kinetic_law(self):

        return self._kinetic_law
    
    def get_list_of_products(self):

        return self._products
    
    def get_list_of_reactants(self):

        return self._reactants