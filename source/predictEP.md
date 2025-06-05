# Analysis EP via RgX

The highly enriched EP links are defined as the following steps:
1. calculate RgX value based on multiomics data
2. Adjustment for RgX by raw value, gene TPM, and chromatin structure
3. Normalized RgX
4. Select maximum distance between gene and sites; 
5. normRgX >1 and RgX ratio > 0.01 as enriched E-P links

## Identify the regulatory sites of a gene set input
selected genes from the E-P pair above, by directly search the same gene ID


## Identify the target gens of regulatory sites 
selected targeted genes from the E-P pairs above, by intersectBed the E-P pairs above

## in silico deletions of given sites

To identify the target genes affected by a given set of genomic loci (especially when they show limited overlap with RgX), we performed an in silico deletion analysis. Specifically, we removed the input loci and recalculated the Rg value for each gene. To assess the statistical significance of these changes, we conducted 100 random deletions of the same number of loci and obtained a background distribution of Rg values. For each gene, we calculated a p-value by comparing the in silico-deleted Rg value with the distribution from the random deletions. The resulting p-values were adjusted for multiple testing using the false discovery rate (FDR), and genes with an FDR-adjusted q-value less than 0.05 were considered significant targets of the input loci.

``` python
import sys
sys.path.append('/home/wang/github/Tichr-CLI/tichr')
from tichr import *
from insilico import *
```

prepare input files
``` python
Rgdir='/home/wang/Tichr/2024Oct-summary/ContextSpecific/RPE_siNIPBL_denovo/resultdf_all_hic/'
RgDF_Ctrl = pd.read_csv(Rgdir+"JQ1minus_BRD4_rep0_RP_RgDf.tsv",header=None,sep="\t")
RgxDF_Ctrl = pd.read_csv(Rgdir+"JQ1minus_BRD4_rep0_RP_RgxDf.tsv",header=None,sep="\t")

deletepeakfile = "/home/wang/Tichr/2025March/predictDEG-withadj/JQ1peak/BRD4.jq1lost.peak"
deletepeakDF=pd.read_csv(deletepeakfile,sep="\t",header=None)
```

- RgDF_Ctrl: Rg file 
- RgxDF_Ctrl: RgX file
- deletepeakfile: insilico deletion file

To validate the result, we can test the lost BRD4 peaks after JQ1 treatment and test if the affected genes are truly DEGs

``` python
silico(RgxDF_Ctrl,deletepeakDF,RgDF_Ctrl,degtype="down",degfdr=0.01,nrandom=20,)
```

the results shows the goodnees of in silico deletions.

<img src="_static/adjrgx/008.png" style="zoom:80%;" />