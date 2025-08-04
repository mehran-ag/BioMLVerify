class BioMLParameter:

    def __init__(self, ID):

        self._ID: str = ID
        self._value: float = None
        self._annotations: dict[str, list[str]] = {}

    @property
    def ID(self):
        return self._ID
    
    @ID.setter
    def ID(self, ID):
        if isinstance(ID, str):
            self._ID = ID
        else:
            raise ValueError("Input foID must be a string")
        

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, value):
        if isinstance(value, (int, float)):
            self._value = value
        else:
            raise ValueError("Input for value must be a number")

    @property
    def annotations(self):
        return self.annotations
    
    @annotations.setter
    def annotations(self, annots):
        if isinstance(annots, dict):
            self._annotations = annots
        else:
            raise ValueError("Input for annotations must be a dictionary")
        


    def get_id(self):

        return self._ID
    
    def get_value(self):

        return self._value