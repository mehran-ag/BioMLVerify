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
    def index(self):
        return self._index
    
    @index.setter
    def index(self, value):
        raise ValueError("Index is read-only")



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
        if isinstance(kfrc, str):
            self._kinetic_forward_rate_constant = kfrc
        else:
            raise ValueError(f"Input for forward kinetic rate constant must be a string!\nEntered:\"{kfrc}\"")
        


    @property
    def kinetic_reverse_rate_constant(self):
        return self._kinetic_reverse_rate_constant
    
    @kinetic_reverse_rate_constant.setter
    def kinetic_reverse_rate_constant(self, krrc):
        if isinstance(krrc, str):
            self._kinetic_reverse_rate_constant = krrc
        else:
            raise ValueError("Input for reverse kinetic rate constant must be a string")
        


    @property
    def kinetic_forward_rate_constant_value(self):
        return self._kinetic_forward_rate_constant_value
    
    @kinetic_forward_rate_constant_value.setter
    def kinetic_forward_rate_constant_value(self, kfrc):
        if isinstance(kfrc, (int, float)):
            self._kinetic_forward_rate_constant_value = kfrc
        else:
            raise ValueError(f"Input for forward kinetic rate constant value must be a number!\nEntered: \"{kfrc}\"")



    @property
    def kinetic_reverse_rate_constant_value(self):
        return self._kinetic_reverse_rate_constant_value
    
    @kinetic_reverse_rate_constant_value.setter
    def kinetic_reverse_rate_constant_value(self, krrc):
        if isinstance(krrc, (int, float)):
            self._kinetic_reverse_rate_constant_value = krrc
        else:
            raise ValueError("Input for reverse kinetic rate constant value must be a number")
      


    @property
    def thermo_forward_rate_constant(self):
        return self._thermo_forward_rate_constant
    
    @thermo_forward_rate_constant.setter
    def thermo_forward_rate_constant(self, tfrc):
        if isinstance(tfrc, str):
            self._thermo_forward_rate_constant = tfrc
        else:
            raise ValueError("Input for thermodynamic rate constant must be a string")



    @property
    def thermo_reverse_rate_constant(self):
        return self._thermo_reverse_rate_constant
    
    @thermo_reverse_rate_constant.setter
    def thermo_reverse_rate_constant(self, trrc):
        if isinstance(trrc, str):
            self._thermo_reverse_rate_constant = trrc
        else:
            raise ValueError("Input for reverse thermodynamic rate constant must be a string")
        


    @property
    def thermo_forward_rate_constant_value(self):
        return self._thermo_forward_rate_constant_value
    
    @thermo_forward_rate_constant_value.setter
    def thermo_forward_rate_constant_value(self, tfrc):
        if isinstance(tfrc, (int, float)):
            self._thermo_forward_rate_constant_value = tfrc
        else:
            raise ValueError("Input for forward thermodynamic rate constant value must be a number")
        


    @property
    def thermo_reverse_rate_constant_value(self):
        return self._thermo_reverse_rate_constant_value
    
    @thermo_reverse_rate_constant_value.setter
    def thermo_reverse_rate_constant_value(self, trrc):
        if isinstance(trrc, (int, float)):
            self._thermo_reverse_rate_constant_value = trrc
        else:
            raise ValueError("Input for reverse thermodynamic rate constant value must be a number")
        


    @property
    def kappa(self):
        return self._kappa
    
    @kappa.setter
    def kappa(self, kp):
        if isinstance(kp, (int, float)):
            self._kappa = kp
        else:
            raise ValueError("Input for Kappa must be a numnber")
        


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
            raise ValueError("Input for kinetic law type must be a string")
        

    @property
    def expanded_kinetic_law(self):
        return self._expanded_kinetic_law
    
    @expanded_kinetic_law.setter
    def expanded_kinetic_law(self, ex_k_law):
        if isinstance(ex_k_law, str):
            self._expanded_kinetic_law = ex_k_law
        else:
            raise ValueError("Input for expanded kinetic law type must be a string")
        


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
        


    @property
    def boundary_condition(self):
        return self._boundary_condition
    
    @boundary_condition.setter
    def boundary_condition(self, BC):
        if isinstance(BC, bool):
            self._boundary_condition = BC
        else:
            raise ValueError("Input for boundary condition must be a boolean")
        


    @property
    def local_parameters(self):
        return self._local_parameters
    
    @local_parameters.setter
    def local_parameters(self, lp):
        if isinstance(lp, list):
            self._local_parameters = lp
        else:
            raise ValueError("Input for local parameters must be a list!")
        
    
    @property
    def klaw_variables(self):
        return self._klaw_variables
    
    @klaw_variables.setter
    def klaw_variables(self, vars):
        if isinstance(vars, list):
            self._klaw_variables = vars
        else:
            raise ValueError("Input for variables must be a list!")
        

    @property
    def mass_action(self):
        return self._mass_action
    
    @mass_action.setter
    def mass_action(self, flag):
        if isinstance(flag, bool):
            self._mass_action = flag
        else:
            raise ValueError("Input for \"mass\" actions must be a Boolean!")