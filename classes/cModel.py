class Model:

    def __init__(self, ID):

        self.ID = ID
        self.reactions = None
        self.species = None
        self.parameters = None

    def getId(self):

        return self.ID

    def getListOfSpecies(self):

        return self.species
    
    def getListOfReactions(self):

        return self.reactions
    
    def getListOfParameters(self):

        return self.parameters