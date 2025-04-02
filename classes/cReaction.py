class Reaction:

    def __init__(self, ID):

        self.ID = ID
        self.annotations = None
        self.reversible: bool = None
        self.kinetic_forward_rate_constant = None
        self.kinetic_reverse_rate_constant = None
        self.thermo_forward_rate_constant = None
        self.thermo_reverse_rate_constant = None
        self.kinetic_law = None
        self.sp_kinetic_law = None
        self.kinetic_law_type = None
        self.reactants = None
        self.products = None

    def getId(self):

        return self.ID
    
    def getReversible(self):

        return self.reversible
    
    def getKineticLaw(self):

        return self.kinetic_law
    
    def getListOfProducts(self):

        return self.products
    
    def getListOfReactants(self):

        return self.reactants