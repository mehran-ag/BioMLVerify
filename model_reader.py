import libsbml
import sys
import os



class BioModel(object):
    '''
    This class takes the file name as a path to read it
    '''

    def __init__(self):
        '''
        Initializes the ModelReader by reading path for a file or a container of files
        '''
        
        self.file_path = None
        self.file_name= None
        self.file_format = None


    def read_file(self, file_path):
                
        
        try:
            if not isinstance(file_path, str):
                raise TypeError
            
            if not os.path.exists(file_path):
                raise FileNotFoundError
            
            self.file_path = file_path
            self.file_name = os.path.basename(file_path)
            self.file_format = os.path.splitext(file_path)[1][1:]
            
        except TypeError:
            print("File not read: Invlaid file type")
            print("\nPlease use a string to enter file path.")

        except FileNotFoundError:
            print("File not read")
            print(f"\nFile not found: {file_path}")
        
        except Exception as e:
            print("File not read")
            print(f"Unexpected error: {e}")

    def _SBML_reader(self):

        reader = libsbml.SBMLReader()
        document = reader.readSBML(self.file_path)
        if document.getNumErrors() > 0:
            print(f"\nThere is (are) {document.getNumErrors()} error(s) in the SBML file.")
            print("\nModel not read")
        else:
            model = document.getModel()
