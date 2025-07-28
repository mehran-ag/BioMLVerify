import libsbml


class BioMLFunctionDefinition():

    def __init__(self, sbml_function_definition):

        if not isinstance(sbml_function_definition, libsbml.FunctionDefinition):
            raise ValueError("Wrong initialization value input for \"FunctionDefinition Class\"")

        self._sbml_function_definition = sbml_function_definition

        self._name = self._sbml_function_definition.getName()

        self._ID = self._sbml_function_definition.getId()

        self._formula = libsbml.formulaToL3String(self._sbml_function_definition.getBody())

        self._arguments = [self._sbml_function_definition.getArgument(n).getName()
                            for n in range(self._sbml_function_definition.getNumArguments())]
        


    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, nm):
        if isinstance(nm, str):
            self._name = nm
        else:
            raise ValueError("Input for name of a FunctionDefinition instance must be string!")
        
    @property
    def ID(self):
        return self._ID
    
    @ID.setter
    def ID(self, i):
        if isinstance(i, str):
            self._ID = i
        else:
            raise ValueError("Input for ID of a FunctionDefinition instance must be string!")
        
    @property
    def formula(self):
        return self._formula
    
    @formula.setter
    def formula(self, fm):
        if isinstance(fm, str):
            self._formula = fm
        else:
            raise ValueError("Input for formula of a FunctionDefinition instance must be string!")
        
    @property
    def arguments(self):
        return self._arguments
    
    @arguments.setter
    def arguments(self, args):
        if isinstance(args, list):
            self._arguments = args
        else:
            raise ValueError("Input for arguments of a FunctionDefinition instance must be list!")
