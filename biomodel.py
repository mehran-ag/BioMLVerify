import libsbml
import sys
import os
from colorama import Fore, Back, Style, init
import matrix_constructor as mc
import exception



class BioModel(object):
    '''
    This class creates a model, either CellML or SBML, that will be processed by other classes
    '''

    def __init__(self):
        '''
        Initializes the ModelReader by reading path for a file or a container of files
        '''
        
        self.file_path = None
        self.file_name= None
        self.file_format = None
        self.model = None
        self.matrix_constructor = mc.MatrixConstructor()


    def read_file(self, file_path):
                
        init( autoreset=True )

        try:
            if not isinstance(file_path, str):
                raise TypeError(Fore.RED + "File path must be a string.")
            
            if not os.path.exists(file_path):
                raise FileNotFoundError(Fore.RED + f"File not found: {file_path}")
            
            self.file_path = file_path
            self.file_name = os.path.basename(file_path)
            self.file_format = os.path.splitext(file_path)[1][1:]
            
        except TypeError as e:
            print(Fore.RED + "File not read!")
            print(Fore.RED + f"Error: {e}")
            return

        except FileNotFoundError as e:
            print(Fore.RED + "File not read!")
            print(Fore.RED + f"Error: {e}")
            return
        
        except Exception as e:
            print(Fore.RED + "File not read!")
            print(Fore.RED + f"Unexpected error: {e}")
            return

        else:
            if self.file_format == 'xml':

                self.model =  self._SBML_reader()

                print(Fore.GREEN + "\n\u27A4\u27A4\u27A4The input file is a SBML model\u27A4\u27A4\u27A4")



    def _SBML_reader(self):

        """
        Reads an SBML file using libSBML.

        :return: SBML model if successful, None otherwise.
        """

        reader = libsbml.SBMLReader()
        document = reader.readSBML(self.file_path)
        if document.getNumErrors() > 0:
            print(f"\nError: The SBML file contains {document.getNumErrors()} error(s).")
            print("\nModel not read")
            return None
        else:
            self.sbmodel = document.getModel()
            return self.sbmodel
        
    def getStoichiometricMatrix(self):

        stoichiometic_matrix = self.matrix_constructor.SBML_stoichiomrtic_matrix_constructor(self.model)

        return stoichiometic_matrix
    
    def getStoichiometricColumnNamesIndices(self):

        return self.matrix_constructor.stoichiometric_matrix_column_names()
    
    def getStoichiometricRowNamesIndices(self):

        return self.matrix_constructor.stoichiometric_matrix_row_names()
    
    def getElementInformationInStoichiometricMatrix(self, i,j):

        return self.matrix_constructor.elementinformation(i,j)
    
    def getThermoConversionMatrix(self):

        try:
            self.matrix_constructor.converstion_matrix_constructor(self.model)
            
        except exception.NoModel as e:
            print(Fore.BLUE + "\nAn error has been raised in \"getThermoConversionMatrix\" function")
            print(Fore.RED + f"{e}")

