from classes.ModelPropertiesMixin import *

class Model(ModelPropertiesMixin):

    def __init__(self, ID):

        self._ID = ID
        self._reactions = None
        self._species = None
        self._parameters = None
        self._function_definitions = None
        self._kinetic_rate_constants_vector: np.ndarray = None
        self._reaction_indices: dict = None
        self._species_indices: dict = None

    def getId(self):

        return self._ID

    def getListOfSpecies(self):

        return self._species
    
    def getListOfReactions(self):

        return self._reactions
    
    def getListOfParameters(self):

        return self._parameters
    
    def setID(self, new_ID):
        self._ID = str(new_ID)