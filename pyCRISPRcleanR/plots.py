import rpy2.rinterface
rpy2.rinterface.set_initoptions(('rpy2', '--no-save', '--no-restore', '--quiet'))
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

d = {'package.dependencies': 'package_dot_dependencies',
     'package_dependencies': 'package_uscore_dependencies'}

base = importr('base', robject_translations=d)
dnacopy = importr("DNAcopy", robject_translations=d)
grdevices = importr('grDevices')
graphics = importr('graphics')


class PlotData(object):

    def __init__(self):
        super().__init__()

    @staticmethod
    def plot_segments(cbs_fc, cbs_normfc, sample_id):
        pdf_prm = {'file': "{}_raw_vs_postCrispr_FC.pdf".format(sample_id), 'width': 7.5, 'height': 7.5}
        grdevices.pdf(**pdf_prm)
        r.par(mfrow=r.c(2, 1))
        for chr_name, (_, _, cnseg_raw) in cbs_fc.items():
            (_, _, cnseg_norm) = cbs_normfc[chr_name]
            plot_prm = {'main': "{}_FCs_chr{}".format(sample_id, chr_name), 'xlab': 'sgRNA_Index',
                        'ylab': 'FCs'}
            dnacopy.plotSample(cnseg_raw, **plot_prm)
            # plot normalised fold changes
            plot_prm = {'main': "{}_postCRISPRcleanR_chr{}".format(sample_id, chr_name), 'xlab': 'sgRNA_Index',
                        'ylab': 'FCs'}
            dnacopy.plotSample(cnseg_norm, **plot_prm)
        grdevices.dev_off()

    @staticmethod
    def box_plot_r(df, title='mytitle', saveto='./myfile.pdf', ylabel='ylabel', xlabel='xlabel'):
        pdf_prm = {'file': saveto + '.pdf', 'width': 7.5, 'height': 7.5}
        grdevices.pdf(**pdf_prm)
        rtbl = pandas2ri.py2ri(df)
        plot_prm = {'main': title, 'names': df.columns, 'xlab': xlabel, 'ylab': ylabel}
        graphics.boxplot(rtbl, **plot_prm)
        grdevices.dev_off()
