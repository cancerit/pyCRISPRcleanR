from abc import ABC, abstractmethod


class AbstractCrispr(ABC):
    '''
      abstract class inilializes files to be comapred and implements two required methods
      to check user input and load config file
    '''

    def __init__(self, **kwargs):
        self.countfile = kwargs['countfile']
        self.libfile = kwargs['libfile']
        self.minreads = kwargs.get('minreads', 30)
        self.mingenes = kwargs.get('mingenes', 3)
        self.outdir = kwargs.get('outdir')
        self.ncontrols = kwargs.get('ncontrols', 1)
        self.ignored_genes = kwargs.get('ignored_genes', [])
        self.runcrispr = kwargs.get('crispr_cleanr', None)
        self.num_processors = kwargs.get('num_processors', 1)
        self.run_mageck = kwargs.get('run_mageck', None)
        self.gene_sig_dir = kwargs.get('gene_signatures', None)
        super().__init__()

    @abstractmethod
    def check_input(self):
        pass
