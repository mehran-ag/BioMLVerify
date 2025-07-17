from classes.BioMLModelPropertiesMixin import *

class BioMLModel(BioMLModelPropertiesMixin):

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
        self._is_mass_action: bool = None

    def get_id(self):

        return self._ID
    
    def get_list_of_compartments(self):

        return self._compartments

    def get_list_of_species(self):

        return self._species
    
    def get_list_of_reactions(self):

        return self._reactions
    
    def get_list_of_parameters(self):

        return self._parameters
    
    def get_list_of_function_definitions(self):

        return self._function_definitions
    
    def set_id(self, new_ID):
        self._ID = str(new_ID)