class Parameter:

    def __init__(self, ID):

        self._ID = ID
        self._annotations = None

    @property
    def ID(self):
        return self._ID
    
    @ID.setter
    def ID(self, ID):
        self._ID = ID

    @property
    def annotations(self):
        return self.annotations
    
    @annotations.setter
    def annotations(self, annots):
        if isinstance(annots, list):
            self._annotations = annots
        else:
            raise ValueError("Input for annotations must be a list")

    def getId(self):

        return self._ID