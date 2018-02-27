CBS_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the CBS data table (SEG)

%(rlibpath)s
library('DNAcopy')

write("Loading gRNA foldchanges", stderr())
tbl = read.delim("%(chr_fc)s")
print(head(tbl))
cna = CNA(cbind(tbl$avgFC), tbl$CHR, tbl$BP,
          data.type="logratio", presorted=T)
write("Smoothing the data", stderr())
smoothed_cna <- smooth.CNA(cna)
write("Segmenting the gRNA data", stderr())
set.seed(0xA5EED)
segment_smoothed_cna = segment(smoothed_cna, verbose=1)
write("Printing the CBS table to standard output", stderr())
segment_smoothed_cna$output$ID = '%(sample_id)s'
write.table(segment_smoothed_cna$output, '', sep='\t', row.names=FALSE)
"""
