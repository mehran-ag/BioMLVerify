from sympy import Basic

class ReactionPropertiesMixin:

    @property
    def ID(self):
        '''Getter for ID'''
        return self._ID
    
    @ID.setter
    def ID(self, ID):
        '''Setter for ID'''
        self._ID = str(ID)

    @property
    def annotations(self):
        '''Getter for annotations'''
        return self._annotations
    
    @annotations.setter
    def annotations(self, annotations):
        '''Setter for annotation - Ensures it is a list'''
        if not isinstance(annotations, list):
            raise ValueError("The annotations must be stored in a list")
        self._annotations = annotations

    @property
    def reversible(self):
        '''Getter for reversible'''
        return self._reversible
    
    @reversible.setter
    def reversible(self, revers):
        if isinstance(revers, bool):
            self._reversible = revers
        else:
            raise ValueError("Input for reversible must be a boolean")
        
    @property
    def kinetic_forward_rate_constant(self):
        return self._kinetic_forward_rate_constant
    
    @kinetic_forward_rate_constant.setter
    def kinetic_forward_rate_constant(self, kfrc):
        if isinstance(kfrc, (int, float)):
            self._kinetic_forward_rate_constant = kfrc
        else:
            raise ValueError("Input for kinetic rate constant must be a number")
        
    @property
    def kinetic_reverse_rate_constant(self):
        return self._kinetic_reverse_rate_constant
    
    @kinetic_reverse_rate_constant.setter
    def kinetic_reverse_rate_constant(self, krrc):
        if isinstance(krrc, (int, float)):
            self._kinetic_reverse_rate_constant = krrc
        else:
            raise ValueError("Input for kinetic rate constant must be a number")
        
    @property
    def thermo_forward_rate_constant(self):
        return self._thermo_forward_rate_constant
    
    @thermo_forward_rate_constant.setter
    def thermo_forward_rate_constant(self, tfrc):
        if isinstance(tfrc, (int, float)):
            self._thermo_forward_rate_constant = tfrc
        else:
            raise ValueError("Input for thermodynamic rate constant must be a number")
        
    @property
    def thermo_reverse_rate_constant(self):
        return self._thermo_reverse_rate_constant
    
    @thermo_reverse_rate_constant.setter
    def thermo_reverse_rate_constant(self, trrc):
        if isinstance(trrc, (int, float)):
            self._thermo_reverse_rate_constant = trrc
        else:
            raise ValueError("Input for thermodynamic rate constant must be a number")
        
    @property
    def kinetic_law(self):
        return self._kinetic_law
    
    @kinetic_law.setter
    def kinetic_law(self, kl):
        if isinstance(kl, str):
            self._kinetic_law = kl
        else:
            raise ValueError("Input for kinetic law must be a string")
        
    @property
    def sp_kinetic_law(self):
        return self._sp_kinetic_law
    
    @sp_kinetic_law.setter
    def sp_kinetic_law(self, skl):
        if isinstance(skl, Basic):
            self._sp_kinetic_law = skl
        else:
            raise ValueError("Input for Sympy Kinetic Law(sp_kinetic_law) must be a Sympy Expression")
        
    @property
    def kinetic_law_type(self):
        return self._kinetic_law_type
    
    @kinetic_law_type.setter
    def kinetic_law_type(self, klt):
        if isinstance(klt, str):
            self._kinetic_law_type = klt
        else:
            raise ValueError("input for kinetic law type must be a string")
        
    @property
    def reactants(self):
        return self._reactants
    
    @reactants.setter
    def reactants(self, reactants):
        if isinstance(reactants, list):
            self._reactants = reactants
        else:
            raise ValueError("Input for reactants must be a list")
        
    @property
    def products(self):
        return self._products
    
    @products.setter
    def products(self, products):
        if isinstance(products, list):
            self._products = products
        else:
            raise ValueError("Input for products must be a list")