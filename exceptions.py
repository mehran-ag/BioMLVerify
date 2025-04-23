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
    Thierror will be raised when initialization of a class fails
    '''
    pass