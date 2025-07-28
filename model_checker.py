import libsbml
import exceptions
import constants as cn
from typing import Union
import sympy as sp
from sympy import symbols
from classes.cBioMLModel import BioMLModel




class ModelChecker(object):
    """
        This class checks the model to see if has the requirements to be eligible for verification checks
        One of the requirements for models is having Mass Action Kinetics
    """


    # ********************************
    # *           Function           *
    # ********************************
    def check_model_reversibility(self, biomlmodel: BioMLModel, return_irreversibles: bool = False) -> Union[bool, tuple[bool, list[str]]]:
        """
            Checks whether all reactions in the input model are reversible.

            If any irreversible reaction is found in the model, the function returns False.
            Optionally, a list of irreversible reactions can also be returned.

            Args:
                biomlmodel (BioMLModel): A model of the BioML class containing species and reactions.
                return_irreversibles (bool): If True, also returns a list of irreversible reactions.

            Returns:
                bool: True if all reactions in the model are reversible, False otherwise.
                list (optional): A list of the IDs of irreversible reactions (returned only if return_irreversibles is True).
        """

        if biomlmodel == None:
            raise exceptions.NoModel("No BioModel has been read!!!")

        biomlmodel_reactions = biomlmodel.get_list_of_reactions()

        irreversible_reactions = []

        reversible = True

        for biomlmodel_reaction in biomlmodel_reactions:

            if biomlmodel_reaction.reversible is not None:

                if not biomlmodel_reaction.reversible:
                    reversible = False
                    irreversible_reactions.append(biomlmodel_reaction.ID)

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
    def _find_variables_in_klaw( self, biomlmodel: BioMLModel ) -> bool:
        """
            Searches for all CellML variables used in all equations in the input BioMLModel and updates the BioMLModel

            Args:
                biomlmodel (BioMLModel): A model of the BioML class containing species and reactions.

            Returns:
                None
        """

        done = True

        biomlmodel_reactions_list = biomlmodel.get_list_of_reactions()

        if biomlmodel_reactions_list:

            for biomlmodel_reaction in biomlmodel_reactions_list:

                if biomlmodel_reaction.expanded_kinetic_law:

                    expanded_kinetic_law_string = biomlmodel_reaction.expanded_kinetic_law

                    variables = ModelChecker._get_variables(expanded_kinetic_law_string)

                    biomlmodel_reaction.klaw_variables = list(dict.fromkeys(variables))

                elif biomlmodel_reaction.kinetic_law:

                    kinetic_law_string = biomlmodel_reaction.kinetic_law

                    variables = ModelChecker._get_variables(kinetic_law_string)

                    biomlmodel_reaction.klaw_variables = list(dict.fromkeys(variables))

                else:

                    done = False

        else:

            done = False
            

        return done

    




    # ********************************
    # *           Function           *
    # ********************************
    def check_mass_action_kinetics(self, biomlmodel: BioMLModel, immediate_return: bool = False) -> bool:
        """
            Checks whether all reactions in the model follow Mass Action Kinetics.

            Each reaction's rate equation is analyzed for conformity with the Mass Action Kinetics pattern.
            Returns True only if all reactions meet this criterion.

            Args:
                biomlmodel (BioMLModel): A model of the BioML class containing species and reactions.
                immediate_return (bool): If True, returns False immediately upon encountering a non-Mass Action reaction.

            Returns:
                bool: True if all reactions follow Mass Action Kinetics, False otherwise.
        """

        flag = True

        if self._find_variables_in_klaw( biomlmodel ):

            biomlmodel_reactions_list = biomlmodel.get_list_of_reactions()

            for biomlmodel_reaction in biomlmodel_reactions_list:

                args = self._make_checking_args( biomlmodel, biomlmodel_reaction )

                status = self._check_kinetic_law(**args)

                if not status:

                    flag = status

                    if immediate_return:
                        return flag

                biomlmodel_reaction.mass_action = status

        else:

            flag = False

        return flag

            








    # ********************************
    # *           Function           *
    # ********************************
    def _make_checking_args(self, biomlmodel: BioMLModel, bioml_reaction: object) -> dict:
        """
            Makes the lists required for checking the Mass Action Kinetics in a reaction.

            Finds all variables used in the kinetic law and stores all ina list.
            Classifies variables used in the kinetic law of the reaction as species, parameters, and compartments.
            Simplifies the kinetic law using Sympy's Simplify function

            Args:
                biomlmodel (BioMLModel): A model of the BioML class containing species and reactions.
                bioml_reaction (BioMLReaction): A reaction intance of BioMLReaction class.

            Returns:
                dict: A dictionary containing:
                    - "species_in_kinetic_law": list of species used in the kinetic law,
                    - "parameters_in_kinetic_law": list of parameters used in the kinetic law,
                    - "compartments_in_kinetic_law": list of compartments used in the kinetic law,
                    - "others_in_kinetic_law": list of other variables not classified as species, parameters, or compartments,
                    - "reactants": list of reactant species,
                    - "products": list of product species,
                    - "kinetic_formula": the original kinetic law expression,
                    - "simp_kinetic_formula": the simplified kinetic law expression,
                    - "klaw_variables": all variables used in the kinetic law.
        """


        reactants = []

        products = []

        species = []

        parameters = []

        compartments = []
        
        for bm_species in biomlmodel.get_list_of_species():

            species.append(bm_species.get_id())

        for bm_parameter in biomlmodel.get_list_of_parameters():

            parameters.append(bm_parameter.get_id())

        for lcl_bm_parameter in bioml_reaction.local_parameters:

            parameters.append(lcl_bm_parameter.get_id())

        for compartment in biomlmodel.get_list_of_compartments():

            compartments.append(compartment)

        parameters = list(dict.fromkeys(parameters))

        klaw_variables = bioml_reaction.klaw_variables

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

        
        for reactant_class in bioml_reaction.get_list_of_reactants():

            reactants.append(reactant_class.get_id())

        for product_class in bioml_reaction.get_list_of_products():

            products.append(product_class.get_id())


        if bioml_reaction.expanded_kinetic_law:

            kinetic_formula = bioml_reaction.expanded_kinetic_law

        elif bioml_reaction.kinetic_law:

            kinetic_formula = bioml_reaction.kinetic_law

        else:

            raise ValueError(f"There is not a kinetic formula for reaction {bioml_reaction.get_id()}")

        # Prepare a space-separated and comma-separated version
        symbol_dict = {klaw_variable: sp.symbols(klaw_variable) for klaw_variable in klaw_variables}

        try:

            symp_kinetic_formula = sp.sympify(kinetic_formula.replace("^", "**"), locals = symbol_dict)

            simp_kinetic_formula = str(sp.simplify(symp_kinetic_formula))

        except:

            simp_kinetic_formula = ''

        return {
            "species_in_kinetic_law": species_in_kinetic_law,
            "parameters_in_kinetic_law": parameters_in_kinetic_law,
            "compartments_in_kinetic_law": compartments_in_kinetic_law,
            "others_in_kinetic_law": others_in_kinetic_law,
            "reactants": reactants,
            "products": products,
            "kinetic_formula": kinetic_formula,
            "simp_kinetic_formula": simp_kinetic_formula,
            "klaw_variables": bioml_reaction.klaw_variables
        }






    # ********************************
    # *           Function           *
    # ********************************
    def _num_klaw_species(self, species_in_kinetic_law: list) -> int:
        """
            Returns the number of species in the kinetic law

            Args:
                species_in_kinetic_law (list): list containing all species of the kinetic law

            Returns:
                int: the number of species in the kinetic law
        """

        return len(species_in_kinetic_law)


            

    # ********************************
    # *           Function           *
    # ********************************
    def _single_product(self, kinetic_law: str, simple_kinetic_law: str) -> bool:
        """
            Checks the input kinetic law to see if it is a single product of terms
            
            Keyword Args:
                kinetic_law (str): kinetic law of a reaction
                simple_kinetic_law (str): a simplified (using sympy simplification command) kinetic law of a reaction
            
            Returns:
               bool: True if is single product of terms, False otherwise.
        """

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
    def _diff_of_products(self, kinetic_law: str, simple_kinetic_law: str) -> bool:
        """
            Checks the input kinetic law to see if it is difference of product of two terms
            
            Keyword Args:
                kinetic_law (str): kinetic law of a reaction
                simple_kinetic_law (str): a simplified (using sympy simplification command) kinetic law of a reaction
            
            Returns:
               bool: True if is difference of product of two terms, False otherwise.
        """

        flag = False

        if self._single_product(kinetic_law, simple_kinetic_law) == False:
            terms = kinetic_law.split("-")
            if len(terms) == 2:
                flag = True

        return flag
    



    # ********************************
    # *           Function           *
    # ********************************
    def _find_frac_parts( self, simp_kinetic_formula: str, klaw_variables: list[str] ) -> dict:
        """
            Finds the numerator and denominator of a fraction in a kinetic law string.

            Keyword Args:
                simp_cellml_eq (str): Simplified CellMl equation using Sympy simplification method
                cellml_vars (list[str]): A list containing CellML variables as strings

            Returns:
                dict: A dictionary with the following keys:
                    - "numerator": A string containing the numerator of the fraction.
                    - "denominator": A string containing the denominator of the fraction.
        """

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
    def _check_kinetic_law(self, **kwargs) -> bool:
        """
            Checks whether the provided kinetic formula follows Mass Action Kinetics.

            Args:
                kwargs (dict): a dictionary containing the following arguments:
                    - "kinetic_formula": A string representing the equation for the kinetic law
                    - "simp_kinetic_formula": A string of the simplified equation of the kinetic law
                    - "species_in_kinetic_law": A list containing the species in the kinetic law
                    - "klaw_variables": A list containing all variables in the kinetic law
                    - "reactants": A list containig all the reactant names of the reaction
                    - "products": A list containig all the product names of the reaction

            Returns:
                bool: True if the equation provided follows the rules of Mass Action equations, False otherwise.
        """


        kinetic_formula = kwargs["kinetic_formula"]
        simp_kinetic_formula = kwargs["simp_kinetic_formula"]
        species_in_kinetic_law = kwargs["species_in_kinetic_law"]
        klaw_variables = kwargs["klaw_variables"]
        reactants = kwargs["reactants"]
        products = kwargs["products"]

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
            fracs = self._find_frac_parts(simp_kinetic_formula, klaw_variables)
            if len(species_in_kinetic_law) > 0:
                for i in range(len(species_in_kinetic_law)):
                    if species_in_kinetic_law[i] in fracs["denominator"]:
                        flag = False

        except:
            pass

        return flag
    









    @staticmethod
    def _identify_variables(ast_node: libsbml.ASTNode, result: list[str]) -> list[str]:
        """
            Recursively traverses an SBML AST (Abstract Syntax Tree) to identify and collect variable names.

            This method inspects the children of a given `ASTNode` and recursively extracts variable names,
            skipping function nodes and continuing traversal for unnamed nodes or functions. It accumulates
            results into the provided `result` list and returns it.

            Parameters:
                ast_node (libsbml.ASTNode): The current AST node to process.
                result (list): A list to accumulate identified variable names during recursion.

            Returns:
                list: A list of variable names (strings) extracted from the AST.

            Raises:
                ValueError: If the recursion depth exceeds `cn.MAX_DEPTH` to prevent infinite recursion.

            Notes:
                - The function relies on a global variable `cur_depth` to track recursion depth.
                - The function uses `getName()` to determine whether a node represents a named variable.
                - Function nodes are further traversed rather than included directly.
                - If `child_node.getName()` returns `None`, it is treated as a non-terminal and is traversed.
        """
    
        global cur_depth
        cur_depth += 1
        
        if cur_depth > cn.MAX_DEPTH:
            raise exceptions.MaxDepth(f"The number of recursions reached maximum depth, {cn.MAX_DEPTH}, but all variables in couldn't be found.")
        
        for idx in range(ast_node.getNumChildren()):
            child_node = ast_node.getChild(idx)

            if child_node.getName() is None:
                additions = ModelChecker._identify_variables(child_node, [])
                result.extend(additions)
            else:
                if child_node.isFunction():
                    additions = ModelChecker._identify_variables(child_node, [])
                    result.extend(additions)
                else:
                    result.append(child_node.getName())
        return result



    @staticmethod
    def _get_variables(kinetic_law_string: str) -> list[str]:
        """
            Extracts all variable names from an SBML ASTNode expression.

            This function initializes recursion depth and invokes a helper method to
            recursively traverse the AST and collect all variable names present in the expression.

            Parameters:
                ast_node (libsbml.ASTNode): The root node of the SBML Abstract Syntax Tree (AST)
                                            representing a mathematical expression.

            Returns:
                list[str]: A list of variable names (strings) extracted from the AST.

            Notes:
                - Uses a global `cur_depth` variable to track recursion depth during traversal.
                - Delegates the recursive extraction to `_identify_variables`.
                - If the root node has a name, it is included in the results.
        """


        global cur_depth
        
        cur_depth = 0

        ast_node = libsbml.parseL3Formula(kinetic_law_string)

        if ast_node is None:
            raise exceptions.NotParsable(f"libsbml.parseL3Formula() couldn't parse the kinetic law, {kinetic_law_string}, for the reaction")
    
        if ast_node.getName() is None:
            variables = []
        else:
            variables = [ast_node.getName()]

        result = ModelChecker._identify_variables(ast_node, variables)

        return result
            