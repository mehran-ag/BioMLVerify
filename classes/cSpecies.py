class Species:

    def __init__(self, ID):
        self.ID = ID
        self.initial_concentration = None
        self.compartment = None
        self.annotations = None
        self.charge = None

    def getId(self):

        return self.ID
    
    def getName(self):

        return self.ID
    
    def getInitialConcentration(self):

        return self.initial_concentration
    
    def getCharge(self):

        return self.charge
    
    def getCompartment(self):

        return self.compartment