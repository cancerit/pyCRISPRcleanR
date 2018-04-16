import plotly.offline as py
import plotly.graph_objs as go
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
    def box_plot_ly(df, title='mytitle', saveto='./myfile.html', ylabel='ylabel', xlabel='xlabel'):
        """
        boxplots for counts data in pandas data frame
        """
        config = {
            'linkText': "Link to documentation plot.ly !!!",
            'scrollZoom': True,
            'displayModeBar': True,
            'editable': True
        }
        layout = go.Layout(title=title, yaxis=dict(title=ylabel), xaxis=dict(title=xlabel))
        data = []
        for col in df.columns:
            data.append(go.Box(y=df[col], name=col, showlegend=False))
        figure = go.Figure(data=data, layout=layout)
        py.plot(figure, filename=saveto + '.html', auto_open=False, config=config)

    @staticmethod
    def plot_segments(cbs_fc, cbs_normfc, sample_id, outdir='./'):
        pdf_prm = {'file': "{}/{}_raw_vs_postCrispr_FC.pdf".format(outdir, sample_id), 'width': 7.5, 'height': 7.5}
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
