import libsbml
import exceptions
import constants as cn
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
        


    def findVariables( self, biomodel ):

        biomodel_reactions_list = biomodel.getListOfReactions()

        for biomodel_reaction in biomodel_reactions_list:

            expanded_kinetic_law_string = biomodel_reaction.expanded_kinetic_law

            variables = ModelChecker.getVariables(expanded_kinetic_law_string)

            biomodel_reaction.variables = list(dict.fromkeys(variables))

            


    









    @staticmethod
    def variableIdentifier(ast_node, result):
    
        global cur_depth
        cur_depth += 1
        
        if cur_depth > cn.MAX_DEPTH:
            raise Exception
        
        for idx in range(ast_node.getNumChildren()):
            child_node = ast_node.getChild(idx)

            if child_node.getName() is None:
                additions = ModelChecker.variableIdentifier(child_node, [])
                result.extend(additions)
            else:
                if child_node.isFunction():
                    additions = ModelChecker.variableIdentifier(child_node, [])
                    result.extend(additions)
                else:
                    result.append(child_node.getName())
        return result



    @staticmethod
    def getVariables(kinetic_law_string):


        global cur_depth
        
        cur_depth = 0

        ast_node = libsbml.parseL3Formula(kinetic_law_string)
    
        if ast_node.getName() is None:
            variables = []
        else:
            variables = [ast_node.getName()]

        result = ModelChecker.variableIdentifier(ast_node, variables)

        return result
            