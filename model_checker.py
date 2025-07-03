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

            if biomodel_reaction.expanded_kinetic_law:

                expanded_kinetic_law_string = biomodel_reaction.expanded_kinetic_law

                variables = ModelChecker.getVariables(expanded_kinetic_law_string)

            else:

                kinetic_law_string = biomodel_reaction.kinetic_law

                variables = ModelChecker.getVariables(kinetic_law_string)


            biomodel_reaction.klaw_variables = list(dict.fromkeys(variables))

    
    def checkMassActionKinetics(self, biomodel):

        self.findVariables( biomodel )

        biomodel_reactions_list = biomodel.getListOfReactions()

        for biomodel_reaction in biomodel_reactions_list:

            self._makeCheckArguments( biomodel, biomodel_reaction )








    def _makeCheckArguments(self, biomodel, bm_reaction):


        reactants = []

        products = []

        species = []

        parameters = []

        compartments = []
        
        for bm_species in biomodel.getListOfSpecies():

            species.append(bm_species.getId())

        for bm_parameter in biomodel.getListOfParameters():

            parameters.append(bm_parameter.getId())

        for lcl_bm_parameter in bm_reaction.local_parameters:

            parameters.append(lcl_bm_parameter.getId())

        for compartment in biomodel.getListOfCompartments():

            compartments.append(compartment)

        parameters = list(dict.fromkeys(parameters))

        klaw_variables = bm_reaction.klaw_variables

        species_in_kinetic_law = []

        parameters_in_kinetic_law = []

        compartments_in_kinetic_law = []

        others_in_kinetic_law = []

        for klaw_variable in klaw_variables:

            if klaw_variable in species:

                species_in_kinetic_law.append(klaw_variable)

            elif klaw_variable in parameters:

                parameters_in_kinetic_law.append(klaw_variable)

            elif klaw_variable in compartments:

                compartments_in_kinetic_law.append(klaw_variable)

            else:

                others_in_kinetic_law.append(klaw_variable)

        
        for reactant_class in bm_reaction.getListOfReactants():

            reactants.append(reactant_class.getId())

        for product_class in bm_reaction.getListOfProducts():

            products.append(product_class.getId())


        if bm_reaction.expanded_kinetic_law:

            kinetic_formula = bm_reaction.expanded_kinetic_law

        elif bm_reaction.kinetic_law:

            kinetic_formula = bm_reaction.kinetic_law

        else:

            raise ValueError(f"There is not a kinetic formula for reaction {bm_reaction.getId()}")
        
        sp_kinetic_formula = kinetic_formula.replace("^", "**")

        

        return {
            "species_in_kinetic_law":species_in_kinetic_law,
            "parameters_in_kinetic_law": parameters_in_kinetic_law,
            "compartments_in_kinetic_law": compartments_in_kinetic_law,
            "others_in_kinetic_law": others_in_kinetic_law,
            "reactants": reactants,
            "products": products,
            "sp_kinetic_formula": sp_kinetic_formula
        }





            


    









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
            