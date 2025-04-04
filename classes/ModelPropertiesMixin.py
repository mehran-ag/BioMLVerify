class ModelPropertiesMixin:

    @property
    def ID(self):
        """Getter for ID"""
        return self._ID

    @ID.setter
    def ID(self, new_ID):
        self._ID = str(new_ID)

    @property
    def reactions(self):
        '''Getter for reactions'''
        return self._reactions
    
    @reactions.setter
    def reactions(self, new_reactions):
        '''Setter for reactions - Ensures it is a list'''
        if not isinstance(new_reactions, list):
            raise ValueError("reactions must be stored in a list")
        self._reactions = new_reactions

    @property
    def species(self):
        '''Getter for species'''
        return self._species
    
    @species.setter
    def species(self, new_species):
        '''Setter for species - Ensures it is a list'''
        if not isinstance(new_species, list):
            raise ValueError("species must be stored in a list")
        self._species = new_species

    @property
    def parameters(self):
        '''Getter for parameters'''
        return self._parameters
    
    @parameters.setter
    def parameters(self, new_parameters):
        '''Setter for parameters - Ensures it is a list'''
        if not isinstance(new_parameters, list):
            raise ValueError("parameters must be stored in a list")
        self._parameters = new_parameters