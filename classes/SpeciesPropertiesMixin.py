class SpeciesPropertiesMixin:

    @property
    def ID(self):
        '''Getter for ID'''
        return self._ID
    
    @ID.setter
    def ID(self, ID):
        '''Setter for ID'''
        self._ID = str(ID)

    @property
    def index(self):
        return self._index
    
    @index.setter
    def index(self, value):
        raise ValueError("Index is read-only")

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        if isinstance(name, str): 
            self._name = name
        else:
            raise ValueError("Name must be a string")

    @property
    def initial_concentration(self):
        return self._initial_concentration
    
    @initial_concentration.setter
    def initial_concentration(self, init_concentration):
        if isinstance(init_concentration, (int, float)):
            self._initial_concentration = init_concentration
        else:
            raise ValueError("Input must be a number")
        
    @property
    def compartment(self):
        return self._compartment
    
    @compartment.setter
    def compartment(self, compartment):
        if isinstance(compartment, str):
            self._compartment = compartment
        else:
            raise ValueError("Input for compartment must be a string!")
        
    @property
    def annotations(self):
        return self._annotations
    
    @annotations.setter
    def annotations(self, annots):
        if isinstance(annots, list):
            self._annotations = annots
        else:
            raise ValueError("Annotations must be stored in a list")
        
    @property
    def charge(self):
        return self._charge
    
    @charge.setter
    def charge(self, charge):
        if isinstance(charge, (int, float)):
            self._charge = charge
        else:
            raise ValueError("Input for Charge must be a number")
        
    @property
    def thermodynamic_rate_constant(self):
        return self._thermodynamic_rate_constant
    
    @thermodynamic_rate_constant.setter
    def thermodynamic_rate_constant(self, trc):
        if isinstance(trc, (int, float)):
            self._thermodynamic_rate_constant = trc
        else:
            raise ValueError("Input must be a number")
        
    @property
    def compound(self):
        return self._compound
    
    @compound.setter
    def compound(self, new_compound):
        if isinstance(new_compound, str):
            self._compound = new_compound
        else:
            raise ValueError("Input must be a string")

    @property
    def composition(self):
        return self._composition
    
    @composition.setter
    def composition(self, new_comp):
        if isinstance(new_comp, dict):
            self._composition = new_comp
        else:
            raise ValueError("Input must be a dictionary mapping elements to their corresponding quantity in the species")
        
    @property
    def chebi_code(self):
        return self._chebi_code
    
    @chebi_code.setter
    def chebi_code(self, new_code):
        if isinstance(new_code, str):
            self._chebi_code = new_code
        else:
            raise ValueError("Input must be a string")