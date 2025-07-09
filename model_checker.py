import libsbml
import exceptions
import constants as cn
from typing import Union
import sympy as sp
from sympy import symbols
import utility




class ModelChecker(object):
    '''
    This class checks the model to see if has the requirements to be eligible for verification checks
    One of the requirements for models is having Mass Action Kinetics
    '''


    # ********************************
    # *           Function           *
    # ********************************
    def check_model_reversibility(self, biomodel, return_irreversibles = False) -> Union[bool, tuple[bool, list]]:

        if biomodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        biomodel_reactions = biomodel.getListOfReactions()

        irreversible_reactions = []

        reversible = True

        for biomodel_reaction in biomodel_reactions:

            if biomodel_reaction.reversible is not None:

                if not biomodel_reaction.reversible:
                    reversible = False
                    irreversible_reactions.append(biomodel_reaction.ID)

            else:

                if return_irreversibles:

                    return None, irreversible_reactions
                
                else:

                    return None

        if return_irreversibles:

            return reversible, irreversible_reactions
        
        else:

            return reversible
        






    # ********************************
    # *           Function           *
    # ********************************
    def _find_variables_in_klaw( self, biomodel ):

        done = True

        biomodel_reactions_list = biomodel.getListOfReactions()

        for biomodel_reaction in biomodel_reactions_list:

            if biomodel_reaction.expanded_kinetic_law:

                expanded_kinetic_law_string = biomodel_reaction.expanded_kinetic_law

                variables = ModelChecker._getVariables(expanded_kinetic_law_string)

                biomodel_reaction.klaw_variables = list(dict.fromkeys(variables))

            elif biomodel_reaction.kinetic_law:

                kinetic_law_string = biomodel_reaction.kinetic_law

                variables = ModelChecker._getVariables(kinetic_law_string)

                biomodel_reaction.klaw_variables = list(dict.fromkeys(variables))

            else:

                done = False
            

        return done

    




    # ********************************
    # *           Function           *
    # ********************************
    def check_mass_action_kinetics(self, biomodel):

        self._find_variables_in_klaw( biomodel )

        biomodel_reactions_list = biomodel.getListOfReactions()

        flag = True

        for biomodel_reaction in biomodel_reactions_list:

            args = self._make_checking_args( biomodel, biomodel_reaction )

            status = self._check_kinetic_law(**args)

            if not status:

                flag = status

            biomodel_reaction.mass_action = status

        return flag

            








    # ********************************
    # *           Function           *
    # ********************************
    def _make_checking_args(self, biomodel, bm_reaction):


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

        # Prepare a space-separated and comma-separated version
        symbol_dict = {klaw_variable: sp.symbols(klaw_variable) for klaw_variable in klaw_variables}

        symp_kinetic_formula = sp.sympify(kinetic_formula.replace("^", "**"), locals = symbol_dict)

        simp_kinetic_formula = str(sp.simplify(symp_kinetic_formula))

        return {
            "species_in_kinetic_law": species_in_kinetic_law,
            "parameters_in_kinetic_law": parameters_in_kinetic_law,
            "compartments_in_kinetic_law": compartments_in_kinetic_law,
            "others_in_kinetic_law": others_in_kinetic_law,
            "reactants": reactants,
            "products": products,
            "kinetic_formula": kinetic_formula,
            "simp_kinetic_formula": simp_kinetic_formula,
            "klaw_variables": bm_reaction.klaw_variables
        }






    # ********************************
    # *           Function           *
    # ********************************
    def _num_klaw_species(self, species_in_kinetic_law):

        return len(species_in_kinetic_law)


            

    # ********************************
    # *           Function           *
    # ********************************
    def _single_product(self, kinetic_law, simple_kinetic_law):

        flag = True

        if "+" in kinetic_law or "-" in kinetic_law:
            flag = False
            if "e-" in kinetic_law or "exp(-" in kinetic_law:
                flag = True
        elif "+" in simple_kinetic_law or "-" in simple_kinetic_law:
            flag = False
            if "e-" in simple_kinetic_law or "exp(-" in simple_kinetic_law:
                flag = True

        return flag
    



    # ********************************
    # *           Function           *
    # ********************************
    def _diff_of_products(self, kinetic_law, simple_kinetic_law):

        flag = False

        if self._single_product(kinetic_law, simple_kinetic_law) == False:
            terms = kinetic_law.split("-")
            if len(terms) == 2:
                flag = True

        return flag
    



    # ********************************
    # *           Function           *
    # ********************************
    def _frac_parts( self, simp_kinetic_formula, klaw_variables ):

        fracs = {"numerator": '', "denominator": ''}

        try:
            symbol_dict = {klaw_variable: sp.symbols(klaw_variable) for klaw_variable in klaw_variables}
            kinetics_sympy_eq = sp.sympify(simp_kinetic_formula, locals=symbol_dict)
            numerator, denominator = kinetics_sympy_eq.as_numer_denom()
            fracs["numerator"] = str(numerator)
            fracs["denominator"] = str(denominator)
        except Exception:
            pass

        return fracs





    # ********************************
    # *           Function           *
    # ********************************
    def _check_kinetic_law(self, **args):

        species_in_kinetic_law = args["species_in_kinetic_law"]
        kinetic_formula = args["kinetic_formula"]
        simp_kinetic_formula = args["simp_kinetic_formula"]
        klaw_variables = args["klaw_variables"]
        reactants = args["reactants"]
        products = args["products"]

        species_in_kinetic_law = species_in_kinetic_law + reactants + products

        species_in_kinetic_law = list(dict.fromkeys(species_in_kinetic_law))

    

        flag = False

        if self._single_product( kinetic_formula, simp_kinetic_formula ) or \
            self._diff_of_products( kinetic_formula, simp_kinetic_formula ):

                if self._num_klaw_species( species_in_kinetic_law ) != 0:

                    flag = True

                    if self._num_klaw_species( species_in_kinetic_law ) == 1:

                        if kinetic_formula.count(species_in_kinetic_law[0]) != 1 or simp_kinetic_formula.count(species_in_kinetic_law[0]) != 1:
                            flag = False
        
        try:
            fracs = self._frac_parts(simp_kinetic_formula, klaw_variables)
            if len(species_in_kinetic_law) > 0:
                for i in range(len(species_in_kinetic_law)):
                    if species_in_kinetic_law[i] in fracs["denominator"]:
                        flag = False

        except:
            pass

        return flag
    









    @staticmethod
    def _variableIdentifier(ast_node, result):
    
        global cur_depth
        cur_depth += 1
        
        if cur_depth > cn.MAX_DEPTH:
            raise Exception
        
        for idx in range(ast_node.getNumChildren()):
            child_node = ast_node.getChild(idx)

            if child_node.getName() is None:
                additions = ModelChecker._variableIdentifier(child_node, [])
                result.extend(additions)
            else:
                if child_node.isFunction():
                    additions = ModelChecker._variableIdentifier(child_node, [])
                    result.extend(additions)
                else:
                    result.append(child_node.getName())
        return result



    @staticmethod
    def _getVariables(kinetic_law_string):


        global cur_depth
        
        cur_depth = 0

        ast_node = libsbml.parseL3Formula(kinetic_law_string)
    
        if ast_node.getName() is None:
            variables = []
        else:
            variables = [ast_node.getName()]

        result = ModelChecker._variableIdentifier(ast_node, variables)

        return result
            