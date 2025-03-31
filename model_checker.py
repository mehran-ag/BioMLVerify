import libsbml
import SBMLKinetics




class ModelChecker(object):
    '''
    This class checks the model to see if has the requirements to be eligible for verification checks
    One of the requirements for models is having Mass Action Kinetics
    '''

    def __init__(self, file_path):

        self.file_path = file_path

    def SBML_checker(self):

        return self._SBML_checker()
    
    def _SBML_checker(self):
        model = SBMLKinetics.readSBML(self.file_path)