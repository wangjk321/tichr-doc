# Context-specific (negative)

## Prepare Rg and RgX file of two conditions

``` python
import sys
sys.path.append('/home/wang/github/Tichr-CLI/tichr')
from tichr import *
from context import *

rgCtrl_path = "negative/Ctrl_RgDf.tsv"
rgTreat_path = "negative/Treat_RgDf.tsv"
rgxCtrl_path = "negative/Ctrl_RgxDf.tsv"
rgxTreat_path = "negative/Treat_RgxDf.tsv"
```

File requirement for RgDf. You can prepare column1~9 as the input file of tichr `candidateGeneFile`.
- column1: gene chr
- column2: gene star
- column3: gene end
- column4: gene ID
- column5: gene symbol
- column6: gene strand
- column7: gene Fold change
- column8: gene FDR
- column9: gene TPM (average of treat and ctrl)
- column10: gene Rg score

File requirment for RgxDf
- column1: site chr
- column2: site start
- column3: site end
- column4: site epigenome signal
- column5: gene ID
- column6: gene chr
- column7: gene start
- column8: gene end
- column9: gene strand
- column10: gene symbol
- column11: site-to-gene weight
- column12: Rgx score
- column13: Rgx ratio

Combine treat and ctrl

``` python
rg_merged,rgx_merged = mergeDF(rgCtrl_path,rgTreat_path,rgxCtrl_path,rgxTreat_path,minRgx=0.1,minRgxRatio=0.01)
```
`minRgx`: filter the site-to-gene links by RgX value > minRgx
`minRgxRatio`: filter the site-to-gene links by RgX Ratio > minRgxRatio

the obtained rgx_merged is:
- column12: RgX score for Ctrl
- column13: RgX ratio, averaged by Ctrl and Treat
- column14: RgX score for Treat

the obtained rg_merged is:
- column10: Rg score for Ctrl
- column11: Rg score for Treat




## check the consistency between Rg-TPM and ΔRg-logFC(TPM)

``` python
prepare_select_by_rank(rg_merged, basedon="rg",filetype="pandas",
                       title="MCF7 with estrogen treatment",label="of ER-binding",genelabel="(all genes)")
```

- `basedon`: could be "rg" or "deltarg". For example, if rg is selected, the program will selected genes with highly related Rg-TPM, and examine if the changes of them， i.e., ΔRg-logFC(TPM) ， are also related.
- `filetype`: could be "pandas" or "file", for the format of input `rg_merged`. If file, a tab-separated without head is needed.
- `diffrank_cutoff` and `sumrank_cutoff`: cutoff to selected highly correlated genes by diffrank or sumrank
- Other parameters is used for visualizations

There are mainly three part of output
1. scatter plot between Rg and TPM, between ΔRg and logFC(TPM).

<img src="_static/negative/001.png" style="zoom:80%;" />

<img src="_static/negative/002.png" style="zoom:80%;" />

2. Select top correlated genes by Rg-TPM diffrank (left and middle), examine the diffrank of ΔRg-logFC(TPM)

<img src="_static/negative/003.png" style="zoom:80%;" />

3. Select top correlated genes by Rg-TPM sumrank (left and middle), examine the sumrank of ΔRg-logFC(TPM)

<img src="_static/negative/004.png" style="zoom:80%;" />

You can also select by ΔRg-logFC(TPM)

<img src="_static/negative/005.png" style="zoom:80%;" />

<img src="_static/negative/006.png" style="zoom:80%;" />


## Identify negatively functional E-P pairs

Identify negative regulatory site-to-gene links by integrating differential chromatin regulation data (TiChR) and gene expression (DEG) analyses.

1. Compute baseline Pearson correlation between gene-level regulatory signal (Rg) and gene expression (TPM).
2. Select significant site-to-gene links (|logFC| > 0.5, FDR < 0.1) connected to differentially expressed genes (|logFC| > 1, FDR < 0.05).
3. Define candidate negative links where regulatory signal and gene expression changes are opposite (ΔRgX * gene logFC < 0).
4. Modify the regulatory signal (RgX') by inverting ΔRgX and weighting by rank difference between ΔRgX and gene logFC.
5. Recalculate a new gene-level regulatory score (Rg') by aggregating RgX' and check if the correlation with TPM improves.
6. Evaluate if the diffrank (difference in ranks) between Rg' and TPM increases for each gene. If all genes show improvement, stop. Otherwise, retain genes with improved diffrank for the next round.
7. Iterate the process until no further improvement is observed, refining the negative link set at each step.



The rg_merged and rgx_merged are used.

``` python
extractNeg(rg_merged, rgx_merged,showInteration=False,
           corrtype="pearson",filetype="pandas",
           outdir="negative",outname="ERnegative",minRgxRatio=0.01)
```
`filetype`:  could be "pandas" or "file", for the format of input `rg_merged` and `rgx_merged`. If file, a tab-separated without head is needed.
`showInteration`: if show the results of each Interation. Default=False
`corrtype`: could be "spearman" or "pearson". Default="pearson"
`outdir`: output directory 
`outname`: output file prefix
`minRgxRatio`: minimal Rgx ratio for negative s2g links.

There are serverl outout:

1. Text information
``` 
pearson correlation - Before: 0.1094, After: 0.3264, Difference: 0.2170
ERnegative
Number of genes with negative sites: 1108
Percentage of genes  with negative sites: 0.2875681287308591
Number of negative sites: 7012
Percentage of negative sites for degs: 0.22255371822134765
```

2. Tab-separated output files for negtive s2g links. Those files are subset of rgx_merged
    - ERnegative_negSite.tsv
    - ERnegative_negSite_filterRgxRatio.tsv

3. Visualizations

    3.1 Final selected site-to-gene links


<img src="_static/negative/007.png" style="zoom:80%;" />

​	3.2 Final selected site-to-gene links

<img src="_static/negative/008.png" style="zoom:80%;" />

3.3 The adjustment for negative links can significanly improve the relationships between Rg-TPM and ΔRg-logFC(TPM)

<img src="_static/negative/009.png" style="zoom:80%;" />


## If you have replicates for epigenomes
``` python
folder="~/Tichr/2024Oct-summary/ContextSpecific/EPSC_shYY1_denovo/resultdf_shYY1_hic/"
shCtrl_rep1_Rg = folder+"ATAC-seq_EPSC_shCtrl_None_rep1_RP_RgDf.tsv"
shCtrl_rep1_Rgx = folder+"ATAC-seq_EPSC_shCtrl_None_rep1_RP_RgxDf.tsv"
shCtrl_rep2_Rg = folder+"ATAC-seq_EPSC_shCtrl_None_rep2_RP_RgDf.tsv"
shCtrl_rep2_Rgx = folder+"ATAC-seq_EPSC_shCtrl_None_rep2_RP_RgxDf.tsv"
shYy1_rep1_Rg = folder+"ATAC-seq_EPSC_shYy1_None_rep1_RP_RgDf.tsv"
shYy1_rep1_Rgx= folder+"ATAC-seq_EPSC_shYy1_None_rep1_RP_RgxDf.tsv"
shYy1_rep2_Rg = folder+"ATAC-seq_EPSC_shYy1_None_rep2_RP_RgDf.tsv"
shYy1_rep2_Rgx = folder+"ATAC-seq_EPSC_shYy1_None_rep2_RP_RgxDf.tsv"

rgCtrl_file = [shCtrl_rep1_Rg,shCtrl_rep2_Rg,]
rgTreat_file = [shYy1_rep1_Rg,shYy1_rep2_Rg,]
rgxCtrl_file = [shCtrl_rep1_Rgx,shCtrl_rep2_Rgx]
rgxTreat_file = [shYy1_rep1_Rgx,shYy1_rep2_Rgx]

rg_merged,rgx_merged = mergeDFmany(rgCtrl_file,rgTreat_file,rgxCtrl_file,rgxTreat_file,minRgx=0.1,minRgxRatio=0.01)

extractNeg(rg_merged, rgx_merged,showInteration=False,
           corrtype="pearson",filetype="pandas",
           outdir="negative",outname="ERnegative",minRgxRatio=0.01,
           epifdr=True)
           
```

