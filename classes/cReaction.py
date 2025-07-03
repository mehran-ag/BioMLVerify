from classes.ReactionPropertiesMixin import *
from typing import Union


class Reaction(ReactionPropertiesMixin):

    _counter: int = 0 #Class variable to track the index

    def __init__(self, ID):

        self._index = Reaction._counter
        Reaction._counter += 1
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
        self._variables: list = []


    @classmethod
    def getCurrentIndex(cls):
        return cls._counter
    
    @classmethod
    def ResetCounter(cls, new_counter_value):
        cls._counter = new_counter_value  # Reset the class-level counter




    def ResetIndex(self):
        self._index = None


    def AssignIndex(self):
        if self._index != None:
            print(f"This reaction has analready has an index: {self._index}")
        else:
            self._index = Reaction._counter
            Reaction._counter += 1
            print(f"Index \"{self._index}\" has now been assigned to this reaction")




    def getId(self):

        return self._ID
    
    def getKineticForwardRateConstantValue(self):
        
        return self._kinetic_forward_rate_constant_value
    
    def getKineticForwardRateConstant(self):
        
        return self._kinetic_forward_rate_constant

    def getReversible(self):

        return self._reversible
    
    def getKineticLaw(self):

        return self._kinetic_law
    
    def getListOfProducts(self):

        return self._products
    
    def getListOfReactants(self):

        return self._reactants