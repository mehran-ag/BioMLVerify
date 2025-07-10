class EmptyList(Exception):
    '''
    This exception catches the exception raised foran empty list
    '''
    pass



class NoModel(Exception):
    '''
    This exception is raised when there is no model to read for a function
    '''
    pass



class LocalParameterConflict(Exception):
    '''
    This error is raised when local parameters for different reactions are similar
    '''
    pass



class Warning(Exception):
    '''
    This type of error will be considered as a warning and will not stop the execution.
    '''
    pass



class WrongInitializationInput(Exception):
    '''
    This error will be raised when initialization of a class fails
    '''
    pass




class NoReverseRateConstant(Exception):
    '''
    This error will be raised when there is no reverse reaction rate constant in constructing the kinetic constants vector and will be useful when checking thermo compatibility
    '''
    pass

class MaxDepth(Exception):
    '''
    This error catches the error in the recursive functions
    '''
    pass

class NotParsable(Exception):
    '''
    This exception is raised when the string cannot be parsed by libsbml FormulatoL3String function
    '''
    pass