"""Segmentation of copy number values."""
import logging
import pandas as pd
from .. import core
from . import cbs

def do_segmentation(cnarr, save_dataframe=False, rlibpath=None, ignoredGenes=[],min_genes=3):
    """Infer copy number segments from the given coverage table."""
    if not len(cnarr):
        return cnarr
    seg_out = ""
    print(cnarr.shape)
    (output,data,segrows,call)=cbs.runCBS(cnarr)
    processed_seg=process_segment(cnarr,output,data,segrows,call)
    #print(cnseg.names)
    # returns vector of dataframes for "data"    "output"  "segRows" "call"
    #do things
    # newcnarr
    #corrected_cnseg=cbs.runCBS(newcnarr)

    #do things
    return None

def process_segment(cnarr,output,data,segrows,call,ignoredGenes=[],min_genes=3):
    newFC=pd.DataFrame(cnarr.avgFC)
    nGeneInSeg={}
    guides={}
    correction={}
    print(newFC.head())
    #correction=rep(0,len(newFC))
    for segment in segrows.itertuples():
        start=segment.startRow
        end=segment.endRow # pyhton indexing is 0 based add 2 to include last row in the range ??
        idxs=list(range(start,end))
        i=segment.Index
        if i == 121:
            """
            print("index:{},Satr:{} end:{}".format(i,start,end))
            print(cnarr.shape)
            print (idxs)
            print(cnarr.head())
            print(cnarr.genes.iloc[[0]])
            print(cnarr.genes.iloc[[9009]])
            continue
            """

            includedGenes=cnarr.genes.iloc[idxs].unique()
            print(includedGenes)
            includedGuides=(start,end)
            if len(ignoredGenes) > 0:
                nGeneInSeg[i]=len(set(includedGenes - ignoredGenes))
            else:
                nGeneInSeg[i]=len(includedGenes)
            if nGeneInSeg[i] >= min_genes:
                newFC[idxs]= newFC.iloc[idxs] - newFC.iloc[idxs].mean()
                print(newFC[idxs])
                correction[idxs]= newFC.iloc[idxs].apply(ns_of_mean)
            guides[i]=includedGuides
            print(guides)
            print(correction)
    #(regions,data2,segrows2,call2)=cbs.runCBS(newFC)
    #correctedFC=newFC
    #regions['gRNA']=guides
    #colnames(regions)<-c('CHR','startp','endp','n.sgRNAs','avg.logFC','guideIdx')
    #res=[correctedFCs=gwSortedFCs,regions=regions].tuple
    return None


def ns_of_mean(row):
    return -(numpy.sign(row.mean(1)))
