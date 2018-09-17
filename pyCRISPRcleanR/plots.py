from math import log10
import plotly.graph_objs as go
import plotly.offline as py
import pandas as pd
import numpy as np
import plotly.tools as tls
from rpy2.robjects import r, pandas2ri
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import plotly.figure_factory as ff
from scipy import stats
import re
import sys
import os

pandas2ri.activate()
numpy2ri.activate()

d = {'package.dependencies': 'package_dot_dependencies',
     'package_dependencies': 'package_uscore_dependencies'}
base = importr('base', robject_translations=d)
dnacopy = importr("DNAcopy", robject_translations=d)
grdevices = importr('grDevices')
pROC = importr('pROC')
PRROC = importr('PRROC')
graphics = importr('graphics')


class PlotData(object):

    def __init__(self):
        super().__init__()

    @staticmethod
    def plotly_conf(cfprm=None):
        """
        default configuration for plotly
        :return:
        """
        config = {
            'linkText': "Link to edit plot (!!!WARNING your data is copied to plotly server) !!!",
            'scrollZoom': True,
            'editable': True
        }

        if cfprm:
            for key, val in cfprm.items():
                config[key] = val
        return config

    @staticmethod
    def histogram_ly(df, title='mytitle', saveto='./myfile', ylabel='ylabel', xlabel='xlabel'):
        """
        :param df:
        :param title:
        :param saveto:
        :param ylabel:
        :param xlabel:
        :return:
        """
        figure = tls.make_subplots(rows=len(df.columns), cols=1, shared_xaxes=True, shared_yaxes=True)
        count = 1
        for col in df.columns:
            x_val = df[col].tolist()
            trace1 = go.Histogram(
                name=col,
                x=x_val,
                opacity=0.75
            )
            figure.append_trace(trace1, count, 1)
            # figure['layout']['xaxis' + str(count)].update(title='')
            # figure['layout']['yaxis' + str(count)].update(title='')
            count += 1
            # if count==len(df.columns):
            #    figure['layout'][xaxis].update(title='counts')
            #    figure['layout'][yaxix].update(title='bins')

        figure['layout'].update(height=700, width=1200, title=title)
        py.plot(figure, filename=saveto + '.html', auto_open=False, config=PlotData.plotly_conf())
        return None

    @staticmethod
    def correlation_matrix_plot(df, title='mytitle', saveto='./myfile', ylabel='ylabel', xlabel='xlabel'):
        figure = ff.create_scatterplotmatrix(df, diag='histogram', height=800, width=800)
        py.plot(figure, filename=saveto + '.html', auto_open=False, config=PlotData.plotly_conf())
        return None

    @staticmethod
    def correlation_plot_ly(df, title='mytitle', saveto='./myfile', ylabel='ylabel', xlabel='xlabel'):
        """
        :
        :param df:
        :param title:
        :param saveto:
        :param ylabel:
        :param xlabel:
        :return:
        """
        dimensions = []
        for col in df.columns:
            d1 = {'label': col,
                  'values': df[col]}
            dimensions.append(d1)

        trace1 = go.Splom(dimensions=dimensions, diagonal=dict(visible=False))
        trace1['dimensions'][1].update(visible=True)
        trace1['showupperhalf'] = False
        annotation_list = []
        yaxis_val = 1
        xcounter = 0
        for col1 in df.columns:
            x = df[col1]
            xcounter += 1
            ycounter = 0
            for col2 in df.columns:
                y = df[col2]
                ycounter += 1

                if xcounter == ycounter:
                    continue
                if xcounter < ycounter:
                    continue
                slope, intercept, r_value, p_value, std_err = stats.linregress(df[col1], df[col2])
                # line = slope * df[col1] + intercept
                format_r_value = "<i>{}_vs_{}: R<sup>2</sup>={:02.2f}</i>".format(col1, col2, r_value)
                annot_dict = dict(x=1,
                                  y=yaxis_val,
                                  xref='paper',
                                  yref='paper',
                                  text=format_r_value,
                                  showarrow=False,
                                  font=dict(size=10)
                                  )
                annotation_list.append(annot_dict)
                yaxis_val -= 0.05

        layout = go.Layout(
            title=title,
            dragmode='select',
            width=600,
            height=600,
            autosize=False,
            hovermode='closest',
            plot_bgcolor='rgba(240,240,240, 0.95)',
            annotations=annotation_list
        )

        figure = dict(data=[trace1], layout=layout)
        py.plot(figure, filename=saveto + '.html', auto_open=False, config=PlotData.plotly_conf())

        return None

    @staticmethod
    def box_plot_ly(df, title='mytitle', saveto='./myfile', ylabel='ylabel', xlabel='xlabel'):
        """
        :param df:
        :param title:
        :param saveto:
        :param ylabel:
        :param xlabel:
        :return:
        """
        data = []
        for col in df.columns:
            data.append(go.Box(y=df[col], name=col, showlegend=False))

        layout = go.Layout(title=title, yaxis=dict(title=ylabel), xaxis=dict(title=xlabel))
        figure = go.Figure(data=data, layout=layout)
        py.plot(figure, filename=saveto + '.html', auto_open=False, config=PlotData.plotly_conf())

    @staticmethod
    def plot_segments(cbs_fc, cbs_normfc, sample_id, outdir='./'):
        """
        :param cbs_fc:
        :param cbs_normfc:
        :param sample_id:
        :param outdir:
        :return:
        """
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

    @staticmethod
    def roc_curve(df, data_type='', fdrth=0.05, saveto='./roc_curve'):
        """
        :param df: data frame with fold changes
        :param data_type:
        :param fdrth:
        :param saveto:
        :return:
        """
        df, roc_auc, sens, _ = PlotData._roc_curve_r(df.tf.values, df.avgFC.values)

        recall = df.sensitivity.values
        tnr = df.specificity.values

        lw = 2

        trace1 = go.Scatter(x=tnr, y=recall,
                            mode='lines',
                            line=dict(color='darkorange', width=lw),
                            showlegend=False
                            )

        trace2 = go.Scatter(x=[1, 0], y=[0, 1],
                            mode='lines',
                            line=dict(color='black', width=0.5),
                            showlegend=False
                            )

        abline = go.Scatter(x=[0, 1], y=[sens, sens],
                            mode='lines',
                            line=dict(color='black', width=0.5, dash='dash'),
                            showlegend=False
                            )
        labels = {
            'x': [0],
            'y': [0],
            'legendgroup': 'group',  # this can be any string, not just "group"
            'name': "Recall {}%FDR={:04.2f} <br></br> AUC={:04.2f}".format(100 * fdrth, sens, roc_auc),
            'opacity': 0,

        }
        layout = go.Layout(title='Receiver operating characteristic({})'.format(data_type),
                           xaxis=dict(title='TNR',
                                      autorange='reversed',
                                      showline=False,
                                      zeroline=False,
                                      ),
                           yaxis=dict(title='Recall',
                                      showline=True,
                                      ),
                           legend=dict(x=.75,
                                       y=.2
                                       )
                           )

        figure = go.Figure(data=[trace1, trace2, abline, labels], layout=layout)

        py.plot(figure, filename=saveto + '_' + data_type + '.html', auto_open=False,
                config=PlotData.plotly_conf())

        return None

    @staticmethod
    def _roc_curve_r(observations, predictions, FDRth=0.05):
        """
        :param observations: known truth set
        :param predictions: all data
        :param FDRth:
        :return:
        """
        obs_rtbl = numpy2ri.py2ri(observations)
        prd_rtbl = numpy2ri.py2ri(predictions)
        roc_prm = {'direction': '>'}
        RES = pROC.roc(obs_rtbl, prd_rtbl, **roc_prm)
        auc = pandas2ri.ri2py(RES.rx2('auc'))[0]
        columns = ['threshold', 'ppv', 'sensitivity', 'specificity']
        coor_prm = {'ret': r.c('threshold', 'ppv', 'sensitivity', 'specificity')}
        COORS = pROC.coords(RES, 'all', **coor_prm)
        cords = numpy2ri.ri2py(COORS)
        df = pd.DataFrame(cords.T, columns=columns)
        FDR5percTh = (df[df.ppv >= (1 - FDRth)])['threshold'].max()
        index_min = min(df[df.threshold <= FDR5percTh].index.tolist())
        threshold = df.at[index_min, 'threshold']
        SENS = df.at[index_min, 'sensitivity']
        SPEC = df.at[index_min, 'specificity']
        return df, auc, SENS, FDR5percTh

    @staticmethod
    def pr_rc_curve(df, data_type='sgRNA', FDRth=0.05, saveto='./pr_rc_curve'):
        """
        :param df: data frame with fold changes
        :param data_type:
        :param FDRth:
        :param saveto:
        :return:
        """
        observations = df.tf.values
        predictions = -df.avgFC.values
        df, auc, SENS = PlotData._pr_rc_curve_r(observations, predictions, FDRth=FDRth)
        h1abline = 1 - FDRth
        h2abline = sum(observations) / observations.size

        trace1 = go.Scatter(x=df.recall.values, y=df.precision.values,
                            mode='lines',
                            line=dict(color='navy', width=2),
                            showlegend=False
                            )

        ablineh1 = go.Scatter(x=[0, 1], y=[h1abline, h1abline],
                              mode='lines',
                              line=dict(color='black', width=0.5, dash='dash'),
                              showlegend=False
                              )
        ablineh2 = go.Scatter(x=[0, 1], y=[h2abline, h2abline],
                              mode='lines',
                              line=dict(color='black', width=0.5),
                              showlegend=False
                              )

        ablinev1 = go.Scatter(x=[SENS, SENS], y=[0, 1],
                              mode='lines',
                              line=dict(color='black', width=0.5),
                              showlegend=False
                              )

        labels = {
            'x': [0],
            'y': [0],
            'legendgroup': 'group',  # this can be any string, not just "group"
            'name': "Recall {}%FDR={:04.2f} <br></br> AUC={:04.2f}".format(100 * FDRth, SENS, auc),
            'opacity': 0,

        }
        layout = go.Layout(title='Precision-Recall({}): AUC={:.2f}'.format(data_type, auc),
                           xaxis=dict(title='Recall'
                                      ),
                           yaxis=dict(title='Precision'),

                           legend=dict(x=.1,
                                       y=.1
                                       )
                           )

        figure = go.Figure(data=[trace1, ablineh1, ablineh2, ablinev1, labels], layout=layout)

        py.plot(figure, filename=saveto + '_' + data_type + '.html', auto_open=False,
                config=PlotData.plotly_conf())

        return None

    @staticmethod
    def _pr_rc_curve_r(observations, predictions, FDRth=0.05):
        """
        :param observations: known truth set
        :param predictions: all data
        :param FDRth:
        :return:
        """
        obs_rtbl = numpy2ri.py2ri(observations)
        prd_rtbl = numpy2ri.py2ri(predictions)
        curve_prm = {'scores.class0': prd_rtbl, 'weights.class0': obs_rtbl, 'curve': True, 'sorted': True}
        prc = PRROC.pr_curve(**curve_prm)
        auc = prc.rx2('auc.integral')[0]
        curve = numpy2ri.ri2py(prc.rx2('curve'))
        cols = ['recall', 'precision', 'threshold']
        df = pd.DataFrame(curve, columns=cols)
        FDR5percTh = - (df[df.precision >= (1 - FDRth)])['threshold'].min()
        index_min = min(df[df.precision >= (1 - FDRth)].index.tolist())
        SENS = df.at[index_min, 'recall']
        threshold = -FDR5percTh

        return df, auc, SENS

    @staticmethod
    def depletion_profile_with_gene_signature(FCsprofile, signatures, df, fdrth=0.05, data_type='genes',
                                              saveto='./genes_depletion_profile'):
        """
        :param FCsprofile: fold change profile
        :param fdrth:
        :param signatures:
        :param df: original data frame
        :param data_type:
        :param saveto:
        :return:
        """
        data_list = PlotData._depletion_profile_with_gene_signature_r(FCsprofile, signatures,
                                                                      df, sigtype=data_type)
        x_fc = data_list[0]
        y_fc = data_list[1]
        FDRpercRANK = data_list[2]
        sig_index_pos_dict_above_fdr = data_list[3]
        sig_index_pos_dict_below_fdr = data_list[4]
        fdr_percent = data_list[5]
        list_b = data_list[6]

        save_image = None
        annotations_list = []
        trace_list = []
        min_log = int(round(min(x_fc) - .5))
        max_log = int(round(max(x_fc) + 1.5))
        labels = [i for i in range(min_log, max_log, 1)]
        tickvalue = [i for i in range(min_log, max_log, 1)]

        trace1 = go.Scatter(x=x_fc, y=list(range(y_fc + 1)),
                            name="",
                            mode='markers',
                            showlegend=False,
                            text=list_b,
                            marker=dict(opacity=1,
                                        color='black'
                                        )
                            )
        trace_list.append(trace1)
        count = max_log
        for sig_name in sig_index_pos_dict_above_fdr:
            # draw scatter above the FDR
            labels.append("{}%".format(fdr_percent[sig_name], '%'))
            sig_name_no_genes = sig_name.split('_genes', -1)[0]
            txt_label = dict(
                x=count,
                y=-0.2,
                xref='x',
                yref='y',
                text=sig_name_no_genes,
                showarrow=False,
                textangle=-90,
                align='left',
            )

            annotations_list.append(txt_label)

            x_range_above = [count] * len(sig_index_pos_dict_above_fdr[sig_name])
            above_trace = go.Scatter(x=x_range_above, y=[i for i, v in sig_index_pos_dict_above_fdr[sig_name]],
                                     showlegend=False,
                                     mode='markers',
                                     marker=dict(symbol="line-ew-open",
                                                 size=19,
                                                 opacity=1,
                                                 color='gray',
                                                 )
                                     )
            trace_list.append(above_trace)
            # draw scatter below FDR
            x_range_below = [count] * len(sig_index_pos_dict_below_fdr[sig_name])
            below_trace = go.Scatter(x=x_range_below, y=[i for i, v in sig_index_pos_dict_below_fdr[sig_name]],
                                     showlegend=False,
                                     mode='markers',
                                     text=sig_name,
                                     textposition='top center',
                                     marker=dict(symbol="line-ew-open",
                                                 size=19,
                                                 opacity=1,
                                                 color='blue'
                                                 )
                                     )
            trace_list.append(below_trace)
            tickvalue.append(count)
            count += 1

        ablineh1 = go.Scatter(x=[min(x_fc) - 1, count], y=[FDRpercRANK, FDRpercRANK],
                              mode='lines',
                              line=dict(color='red', width=2.5, dash='dash'),
                              showlegend=False,
                              )
        trace_list.append(ablineh1)
        # axis and labels
        ypos = int(round(log10(FDRpercRANK)))
        fdr_label = dict(
            x=max_log - 1, y=ypos + 0.12,
            xref='x',
            yref='y',
            text="{}%FDR".format(100 * fdrth),
            showarrow=False,
            font=dict(
                size=16,
                color='red'
            ),
            align='center'
        )

        x_label1 = dict(
            x=0.25,
            y=-0.075,
            showarrow=False,
            text="LogFC",
            xref='paper',
            yref='paper',
            font=dict(size=14)
        )

        x_label2 = dict(
            x=0.85,
            y=-0.075,
            showarrow=False,
            text="% of genes above and below {}%FDR cutoff".format(100 * fdrth),
            xref='paper',
            yref='paper',
            font=dict(size=14)
        )

        annotations_list.append(fdr_label)
        annotations_list.append(x_label1)
        annotations_list.append(x_label2)

        layout = go.Layout(
            title='Depletion Profile: {}'.format(data_type),
            xaxis=dict(
                title="",
                showline=True,
                showgrid=False,
                zeroline=True,
                range=[min_log, count],
                ticks="",
                showticklabels=True,
                ticktext=labels,
                tickvals=tickvalue,
                dtick=1
            ),
            yaxis=dict(title='Depeltion Rank',
                       autorange='reversed',
                       showgrid=False,
                       range=[1, y_fc],
                       type='log',
                       showline=True,
                       zeroline=False,
                       scaleanchor="x",
                       scaleratio=0.2,
                       domain=[1, 1],
                       dtick=1,
                       ),
            annotations=annotations_list

        )

        figure = go.Figure(data=trace_list, layout=layout)

        if save_image:
            py.plot(figure, filename=saveto + '_' + data_type + '.html', auto_open=False,
                    config=PlotData.plotly_conf(),
                    image='jpeg', image_filename=saveto, image_width=1200, image_height=800)
        else:
            py.plot(figure, filename=saveto + '_' + data_type + '.html', auto_open=False,
                    config=PlotData.plotly_conf())

        return None

    @staticmethod
    def _depletion_profile_with_gene_signature_r(FCsprofile_np, signatures, df, sigtype='genes'):
        """
        :param FCsprofile_np: profile numpy array
        :param signatures: signatures dictionary
        :param df: data frame
        :param sigtype: type of signature data
        :return x_fc, y_fc, FDRpercRANK, sig_index_pos_dict_above_fdr, sig_index_pos_dict_below_fdr,
        fdr_percent, list_b
        """
        print("Draw signature profile")
        FCsprofile = FCsprofile_np.to_frame()
        FCsprofile.sort_values(by=['avgFC'], inplace=True)
        # get ROC data
        df, auc, sens, FDR5percTh = PlotData._roc_curve_r(df.tf.values, df.avgFC.values)
        FDRpercRANK = (FCsprofile[FCsprofile.avgFC <= FDR5percTh]).shape[0]
        # print("{}:{}:{}:FDRpercRANK:{},FDR5percTh:{}".format(FCsprofile.head(), auc,
        # sens,FDRpercRANK,FDR5percTh))

        y_fc = FCsprofile.shape[0]
        x_fc = FCsprofile.avgFC.values

        sig_index_pos_dict_below_fdr = {}
        sig_index_pos_dict_above_fdr = {}
        fdr_percent = {}
        list_b = FCsprofile.index.tolist()
        sigtype_re = re.compile(sigtype + '$')
        for sig in signatures:
            if sigtype_re.search(sig):
                # find positions in the profile
                index_tuple = PlotData._get_indices(list_b, signatures[sig])
                sig_index_pos_dict_below_fdr[sig] = [(i, v) for i, v in index_tuple if i <= FDRpercRANK]
                sig_index_pos_dict_above_fdr[sig] = [(i, v) for i, v in index_tuple if i > FDRpercRANK]
                fdr_percent[sig] = round((len(sig_index_pos_dict_below_fdr[sig]) / len(index_tuple)) * 100)
                # print("Signature:{}:index above:{}: below:{}:fdr_perc:{}".format(sig,
                # len(sig_index_pos_dict_above_fdr[sig]),
                # len(sig_index_pos_dict_below_fdr[sig]),fdr_percent[sig] ))

        return [x_fc, y_fc, FDRpercRANK, sig_index_pos_dict_above_fdr, sig_index_pos_dict_below_fdr,
                fdr_percent, list_b]

    @staticmethod
    def _get_indices(a, b):
        b_set = set(b)
        return [(i, v) for i, v in enumerate(a) if v in b_set]

    @staticmethod
    def impact_on_phenotype(mo_uncorrected_file, mo_corrected_file, fdrth=0.05, saveto='./impact_on_phenotype',
                            exp_name='myexp'):
        """
        :param fdrth: fdr cutoff
        :param saveto: path to save file
        :param exp_name: user defined name for experiment
        :param mo_uncorrected_file: mageck uncorrected output
        :param mo_corrected_file: mageck corrected output
        :return: None
        """
        pre = pd.read_csv(mo_uncorrected_file, compression='infer', sep="\t", index_col='id')
        post = pd.read_csv(mo_corrected_file, compression='infer', sep="\t", index_col='id')
        # print(pre.head())

        pre.sort_index(inplace=True)
        post.sort_index(inplace=True)
        pre_genes = set(pre.index.tolist())
        preD = set(pre[(pre['neg|fdr'] < fdrth) & (pre['pos|fdr'] >= fdrth)].index.tolist())
        preE = set(pre[(pre['neg|fdr'] >= fdrth) & (pre['pos|fdr'] < fdrth)].index.tolist())
        preNULL = pre_genes.difference(preD.union(preE))

        # print("all:{}pred:{}preE{}preNULL:{}".format(len(pre_genes),len(preD),len(preE),len(preNULL)))

        post_genes = set(post.index.tolist())
        postD = set(post[(post['neg|fdr'] < fdrth) & (post['pos|fdr'] >= fdrth)].index.tolist())
        postE = set(post[(post['neg|fdr'] >= fdrth) & (post['pos|fdr'] < fdrth)].index.tolist())
        postNULL = post_genes.difference(postD.union(postE))

        # print("all:{}postd:{}postE{}postNULL:{}".format(len(post_genes), len(postD),
        # len(postE), len(postNULL)))

        aDD = len(preD.intersection(postD))
        aDN = len(preD.intersection(postNULL))
        aDE = len(preD.intersection(postE))

        aND = len(preNULL.intersection(postD))
        aNN = len(preNULL.intersection(postNULL))
        aNE = len(preNULL.intersection(postE))

        aED = len(preE.intersection(postD))
        aEN = len(preE.intersection(postNULL))
        aEE = len(preE.intersection(postE))

        cm = np.matrix([[aDD, aDN, aDE], [aND, aNN, aNE], [aED, aEN, aEE]])
        IMPACTEDg = ((np.triu(cm, 1) + np.tril(cm, -1)) / cm.sum()).sum() * 100
        IMPACTED_phenGenes = (cm[0, 1] + cm[2, 1] + cm[0, 2] + cm[1, 2]) / (cm[[0, 2]]).sum() * 100
        DISTORTEDg = (cm[2, 0] + cm[0, 2]) / cm.sum() * 100
        DISTORTED_phenGenes = (cm[2, 0] + cm[0, 2]) / (cm[[0, 2]]).sum() * 100
        cm_fitmat = np.divide(cm.T, np.concatenate([cm.sum(1).T] * 3)) * 100

        # print("IMPACTEDg:{}IMPACTED_phenGenes:{}DISTORTEDg:{}DISTORTED_phenGenes:{},cm:{}".format(IMPACTEDg,
        # IMPACTED_phenGenes, DISTORTEDg,DISTORTED_phenGenes,cm_fitmat*100))

        original_counts = cm.sum(1).T.tolist()[0]
        label_list = ['loss_of_fitness', 'no_phenotype', 'gain_of_fitness']
        mod_label = ["{}<br>{}</br>".format(original_counts[i], label_list[i]) for i in range(3)]

        barplot_list = []

        for i in range(3):
            trace1 = go.Bar(
                x=mod_label,
                y=cm_fitmat[i, :].tolist()[0],
                name=label_list[i]
            )
            barplot_list.append(trace1)

        layout = go.Layout(
            title=exp_name,

            xaxis=dict(
                title="Original Counts",
            ),
            yaxis=dict(title='% of genes',
                       ),
            barmode='stack',
            legend=dict(
                x=1,
                y=1,
                orientation="v"
            ),
            annotations=[
                dict(
                    x=1.17,
                    y=1.01,
                    align="right",
                    valign="top",
                    text='Corrected Counts',
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    xanchor="top",
                    yanchor="top"
                )]
        )
        figure = go.Figure(data=barplot_list, layout=layout)

        py.plot(figure, filename=saveto + '.html', auto_open=False,
                config=PlotData.plotly_conf(cfprm={'theme': 'ggplot'}))

        pichart_dict = {('Overall impact', 'Rest of the genes'): [IMPACTEDg, 100 - IMPACTEDg],
                        ('Overall distortion', 'Rest of the genes'): [DISTORTEDg, 100 - DISTORTEDg],
                        ('Impact (G/L fitness genes)', 'Rest of the genes'): [IMPACTED_phenGenes,
                                                                              100 - IMPACTED_phenGenes],
                        ('Distortion (G/L fitness genes)', 'Rest of the genes'): [DISTORTED_phenGenes,
                                                                                  100 - DISTORTED_phenGenes]}
        pi_colors = ['red', 'blue', 'green', 'orange']
        pi_x = [[0, .48], [.52, 1], [0, .48], [.52, 1]]
        pi_y = [[0, .49], [0, .49], [.51, 1], [.51, 1]]

        # pie chart
        pie_list = []
        count = 0

        for label_tuple in pichart_dict:
            pie_data = {
                'labels': list(label_tuple),
                'values': pichart_dict[label_tuple],
                'type': 'pie',
                'domain': {'x': pi_x[count],
                           'y': pi_y[count]
                           },
                'marker': {'colors': [pi_colors[count], '#FFFFFF'],
                           'line': {'color': '#000000',
                                    'width': 0.5},
                           },
                'hoverinfo': 'label+percent',
                'textinfo': 'label'
            }
            count += 1

            pie_list.append(pie_data)
        layout = go.Layout(
            title="Impact of Correction on Genes",
        )

        figure = go.Figure(data=pie_list, layout=layout)

        py.plot(figure, filename=saveto + 'piechart.html', auto_open=False,
                config=PlotData.plotly_conf())

        return None
