from abc import ABC, abstractmethod


class AbstractCrispr(ABC):
    '''
      abstract class inilializes files to be comapred and implements two required methods
      to check user input and load config file
    '''

    def __init__(self, **kwargs):
        self.countfile = kwargs['countfile']
        self.libfile = kwargs['libfile']
        self.expname = kwargs.get('expname', 'myexperiment')
        self.minreads = kwargs.get('minreads', 30)
        self.mingenes = kwargs.get('mingenes', 3)
        self.outdir = kwargs.get('outdir', './')
        self.ncontrols = kwargs.get('ncontrols', 1)
        self.sample = kwargs.get('sample', 'mySample')
        self.ignored_genes = kwargs.get('ignored_genes', [])
        self.num_processors = kwargs.get('num_processors', 1)
        super().__init__()

    @abstractmethod
    def check_input(self):
        pass
