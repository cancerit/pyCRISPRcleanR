from abc import ABC, abstractmethod


class AbstractCompare(ABC):
    '''
      abstract class inilializes files to be comapred and implements two required methods
      to check user input and load config file
    '''

    def __init__(self, **kwargs):
        self.file_a = kwargs['archive_a']
        self.file_b = kwargs['archive_b']
        self.json_file = kwargs.get('json_config', None)
        self.cmp_type = kwargs.get('cmp_type', None)
        self.outfile = kwargs.get('outfile', None)
        self.get_config()
        super().__init__()

    @abstractmethod
    def check_input(self):
        pass

    @abstractmethod
    def get_config(self):
        pass
