from classes.ModelPropertiesMixin import *

class Model(ModelPropertiesMixin):

    def __init__(self, ID):

        self._ID = ID
        self._compartments = []
        self._reactions = []
        self._species = []
        self._parameters = []
        self._function_definitions = []
        self._kinetic_rate_constants_vector: np.ndarray = None
        self._reaction_indices: dict = None
        self._species_indices: dict = None

    def getId(self):

        return self._ID
    
    def getListOfCompartments(self):

        return self._compartments

    def getListOfSpecies(self):

        return self._species
    
    def getListOfReactions(self):

        return self._reactions
    
    def getListOfParameters(self):

        return self._parameters
    
    def getListOfFunctionDefinitions(self):

        return self._function_definitions
    
    def setID(self, new_ID):
        self._ID = str(new_ID)