import libsbml
import exceptions
import SBMLKinetics
from typing import Union




class ModelChecker(object):
    '''
    This class checks the model to see if has the requirements to be eligible for verification checks
    One of the requirements for models is having Mass Action Kinetics
    '''


    def check_model_reversibility(self, biomodel, return_irreversibles = False) -> Union[bool, tuple[bool, list]]:

        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        biomodel_reactions = biomodel.getListOfReactions()

        irreversible_reactions = []

        reversible = True

        for biomodel_reaction in biomodel_reactions:

            if not biomodel_reaction.reversible:
                reversible = False
                irreversible_reactions.append(biomodel_reaction.ID)

        if return_irreversibles:

            return reversible, irreversible_reactions
        
        else:

            return reversible
        

        


    def SBML_checker(self):

        return self._SBML_checker()
    
    def _SBML_checker(self):
        model = SBMLKinetics.readSBML(self.file_path)