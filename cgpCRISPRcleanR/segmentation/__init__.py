"""Segmentation of copy number values."""
import logging
import pandas as pd
from .. import core
from . import cbs
import numpy as np
pd.set_option('mode.chained_assignment', 'raise')

def do_segmentation(cnarr, save_dataframe=False, rlibpath=None, ignoredGenes=[],min_genes=3):
    """Infer copy number segments from the given coverage table."""
    if not len(cnarr):
        return cnarr
    seg_out = ""
    print(cnarr.shape)
    (output,data,segrows,calls)=cbs.runCBS(cnarr)
    print("output:{} , data:{} , segrows:{}".format(output.shape,data.shape,segrows.shape))

    (correctedFC,regions)=process_segment(cnarr,output,segrows)
    #print(correctedFC.head())
    #print(correctedFC.tail())
    #print(regions.head())
    #print(regions.tail())
    return (correctedFC,regions)

def process_segment(cnarr,output,segrows,ignoredGenes=[],min_genes=3):
    cnarr.is_copy= False  # to avoid   SettingWithCopyWarning
    #CHR  startp    endp   genes avgFC  BP
    output.is_copy= False
    cnarr['correction']=0
    cnarr['correctedFC']=cnarr.avgFC

    #ID chrom  loc.start  loc.end  num.mark  seg.mean
    output=output.rename(columns={'chrom': 'CHR','loc.start':'startp','loc.end':'endp', 'num.mark':'n.sgRNAs', 'seg.mean':'avg.logFC'})
    output.drop(output.columns[0], axis=1, inplace=True)
    output['startRow']=0
    output['endRow']=0
    output['nGenes']=0
    nGeneInSeg=0
    for segment in segrows.itertuples():
        start=segment.startRow -1 #pyhton indexing is 0 based add 1 to include last row in the range ??
        end=segment.endRow #
        idxs=list(range(start, end))
        i=segment.Index
        #print("index:{},Satr:{} end:{}".format(i,start,end))
        includedGenes=cnarr.genes.iloc[idxs].unique()
        output.iat[i,5]=start
        output.iat[i,6]=end
        if len(ignoredGenes) > 0:
            nGeneInSeg=len(set(includedGenes - ignoredGenes))
        else:
            nGeneInSeg=len(includedGenes)
        output.iat[i,7]=nGeneInSeg
        if nGeneInSeg >= min_genes:
            cnarr.iloc[idxs,cnarr.columns.get_loc('correctedFC')]=cnarr.avgFC.iloc[idxs] - cnarr.avgFC.iloc[idxs].mean()
            cnarr.iloc[idxs,cnarr.columns.get_loc('correction')] = -np.sign(cnarr.correctedFC.iloc[idxs].mean())
    return (cnarr,output)
