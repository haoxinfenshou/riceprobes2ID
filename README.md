# riceprobes2ID

This is an R package for converting rice microarray ID to RAP-ID and merging multi probes for one gene.

Until now, 3 platforms including 3,555 samples are available:

GPL6864: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6864
GPL8852: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8852
GPL2025: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2025

Firstly, an expression matrix is needed, and the rowname is microarray ID and the colname is sample ID(GSMXXXX).
The expression matrix can be downloaded and extracted following NCBI GEO2R's protocols.

Then, rowname with microarray ID can be converted using the function probe_trans().

There are 3 parameters in probe_trans(), 'GPL', 'expr_matrix' and 'merge_by'.

GPL requires a GPL ID, GPL = 'GPL6864', GPL = 'GPL8852' or GPL = 'GPL2025'.

expr_matrix requires an above-mentioned expression matrix.

merge_by can choose a method of merging multi probes for one gene, now only mean value is available.
If there are 2 probes for one gene, the signal intensity of one probe is 2 and another is 4, so after converting and merging, 
we can see in the exported expression matrix and the signal intensity of this gene is 3.

exported_expression_matrix <- riceprobes2ID::probe_trans('GPL8852', expr_matrix = expr)

The exported_expression_matrix can be used for further DEGs calculation.
