from _classes.BioMLModelPropertiesMixin import *

class BioMLModel(BioMLModelPropertiesMixin):

    def __init__(self, ID):

        self._ID: str = ID
        self._compartments: list[str] = []
        self._reactions: list[object] = []
        self._species: list[object] = []
        self._parameters: list[object] = []
        self._function_definitions: list[object] = []
        self._kinetic_rate_constants_vector: np.ndarray = None
        self._reaction_indices: dict = None
        self._species_indices: dict = None
        self._is_mass_action: bool = None
        self._is_direct_conversion: bool = False

        self._element_indices_dict = {}

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

    
    def get_element_indices_dict(self):

        if self._element_indices_dict():

            return self._element_indices_dict()
        
        else:

            try:

                self.mk_element_indices_dict()

                return self._element_indices_dict

            except Exception:
                
                return {}
        


    def mk_element_indices_dict(self):

        counter = 0

        for single_species in self._species:

            if single_species.composition is None:
                raise ValueError(f"Chemical composition is not known for {single_species.name if single_species.name is not None else single_species.ID}\n       Elemental matrix cannot be costructed!")

            for compound in single_species.composition.keys():

                if compound not in self._element_indices_dict:

                    self._element_indices_dict[compound] = counter

                    counter += 1

        return self._element_indices_dict


        


