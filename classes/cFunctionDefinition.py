import libsbml


class FunctionDefinition():

    def __init__(self, sbml_function_definition):

        if not isinstance(sbml_function_definition, libsbml.FunctionDefinition):
            raise ValueError("Wrong initialization value input for \"FunctionDefinition Class\"")

        self.sbml_function_definition = sbml_function_definition

        self.name = self.sbml_function_definition.getName()

        self.ID = self.sbml_function_definition.getId()

        self.formula = libsbml.formulaToL3String(self.sbml_function_definition.getBody())

        self.arguments = [self.sbml_function_definition.getArgument(n).getName()
                            for n in range(self.sbml_function_definition.getNumArguments())]
        


        @property
        def name(self):
            return self.name
        
        @name.setter
        def name(self, nm):
            if isinstance(nm, str):
                self.name = nm
            else:
                raise ValueError("Input for name of a FunctionDefinition instance must be string!")
            
        @property
        def ID(self):
            return self.ID
        
        @ID.setter
        def ID(self, i):
            if isinstance(i, str):
                self.ID = i
            else:
                raise ValueError("Input for ID of a FunctionDefinition instance must be string!")
            
        @property
        def formula(self):
            return self.formula
        
        @formula.setter
        def ID(self, fm):
            if isinstance(fm, str):
                self.formula = fm
            else:
                raise ValueError("Input for formula of a FunctionDefinition instance must be string!")
            
        @property
        def arguments(self):
            return self.arguments
        
        @arguments.setter
        def arguments(self, args):
            if isinstance(args, list):
                self.arguments = args
            else:
                raise ValueError("Input for arguments of a FunctionDefinition instance must be list!")
