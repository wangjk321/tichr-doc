# 2 Command-Line Interface

## Overall design

The core function of **Tichr** -  `compute` and `negative` -  can be used in the form of a command-line tool. 

Unlike the **Tichr** API, the **Tichr** CLI is more convenient to use. It does not require distributed processing and can directly obtain **Rg** and **RgX** results through a single command line. However, in terms of flexibility, it falls short compared to the **Tichr** API. 

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 10px; color: #000000ff;; background-color: #b5cabaff; border-color: #faebcc;">
<p><strong>Important:</strong></p>
<p>
This command-line interface (CLI) provides the most common functions of TICHR.  
If you need more advanced or customized functionalities, please refer to the other sections of the documentation.
</p>
</div>

You can check the framework of tichr by:

``` shell
tichr -h
```

The output will be:

``` text
usage: __main__.py [-h] [-V] {calcu,adjust,neg} ...

TICHR is software to analyse transcriptional regulation by integrating Epigenome
(ChIP-seq etc.), 3D genome (Hi-C) and Transcriptome (RNA-seq). See the project
page at https://github.com/wangjk321/tichr. This command-line provides basic
usages; more functions are available within the Python API.

positional arguments:
  {calcu,adjust,neg}  Choose the mode to use sub-commands
    calcu             Calculate Rg and RgX based on multiomics data
    adjust            Adjust S2G regulation
    neg               Identify context-specific repressive functions

options:
  -h, --help          show this help message and exit
  -V, --version       Show tichr version
```


## 1. Compute RgX and Rg

#### Description 

The usage of **Tichr**  `calculate` command-line tool includes the following parameters:

``` text
$ tichr calcu [-h] [--refGeneFile REFGENEFILE] [--TSSrange TSSRANGE]
                  [--S2Gmax S2GMAX] [--hicfilepath HICFILEPATH]
                  [--readFileList2 READFILELIST2]
                  [--macs2species MACS2SPECIES]
                  [--binResolution BINRESOLUTION]
                  [--blackregion BLACKREGION] [--fixPeakWidth]
                  [--tmpdir TMPDIR] [--coverageMethod COVERAGEMETHOD]
                  [--spmr] [--multiBAMmerge MULTIBAMMERGE]
                  [--file_type FILE_TYPE] [--hicRes HICRES]
                  [--hicDataType HICDATATYPE] [--hicNormType HICNORMTYPE]
                  [--juicertool JUICERTOOL] [--threads THREADS]
                  [--contactNorm CONTACTNORM] [--ifUseHiCRef]
                  [--weightType WEIGHTTYPE]
                  [--fixedFunctionType FIXEDFUNCTIONTYPE]
                  [--halfDistance HALFDISTANCE] [--setpromoter1]
                  [--threadscalcu THREADSCALCU] [--outname OUTNAME]
                  [--outhead] [--noise_ratio NOISE_RATIO]
                  [--noise_quantile NOISE_QUANTILE]
                  readFileList candidateSite candidateGeneFile gtfile
```

The `calcu` method of **Tichr** CLI is consistent with that of the **Tichr** API `compute`. If you are interested in detailed information about certain parameters, you can find more help <a href="https://tichr.readthedocs.io/en/latest/1.Compute.html">here</a>. 

#### Typical scripts 

Here is an example to show you how to apply Tichr `calcu` CLI in a example dataset:

```shell
testdata=./Data/Calculate

candidatesite="$testdata/candidatePeak.bed"
candidateGeneFile="$testdata/tichr_allcoding_hg38.tsv"
DNase_rep1="$testdata/DNase_rep1.bigwig"
DNase_rep2="$testdata/DNase_rep2.bigwig"
gtfile="$testdata/hg38_gt.tsv"

#Use  distance-based weights
tichr calcu $DNase_rep1,$DNase_rep2 \
        $candidatesite $candidateGeneFile $gtfile \
        --S2Gmax 500000 --fixPeakWidth --spmr --file_type bigWig  \
        --weightType fixedFunction --fixedFunctionType exponential\
        --halfDistance 10000 --outname calcuout --outhead

#Use  hic-based weights
hicfilepath="../ENCFF621AIY.hic"
python -m tichr calcu $DNase_rep1,$DNase_rep2 \
        $candidatesite $candidateGeneFile $gtfile \
        --S2Gmax 500000 --hicfilepath $hicfilepath \
        --fixPeakWidth --spmr --file_type bigWig  \
        --hicRes 50000 --hicNormType VC_SQRT --threads 12 \
        --contactNorm 'default' --weightType hic \
        --outname calcuout --outhead
```


#### Output information

The output information is like:

``` text
Creating Tichr object...
***Checking file availablity
Start makeSiteBed...
***Temporary directory created at: tichr_tmp_RRAclGYFOa
***Setting candidate sites
......Use a user-defined candidate bed file
Finish makeSiteBed
Start makeSiteBdg...
Loading genome table: ./Data/Calculate/hg38_gt.tsv
Loading BED file: tichr_tmp_RRAclGYFOa/candidatesite_final.bed
✔ No BED errors found. All intervals are within genome bounds.
***Extracting read counts from epigenome files
use spmr
......Input epigenome signal is Bigwig.
......Finished.
Finish makeSiteBdg
Start process HiC...
***Processing hic ...
......Using 12 threads to process Hi-C. More threads require larger memory.
......Normalizing the contact weight.
......Finished the Hi-C weight process.
Finish process HiC
Start Computing...
***Computing RgX and Rg
***Using hic mode
***Processing genes: 100%|█████████████████| 20049/20049 [00:31<00:00, 632.53gene/s]
***Finished
Saving files...
***Save RgX and Rg to tsv tables.
......Saved to calcuout_RgxDf.tsv.gz and calcuout_RgDf.tsv.gz
Finish Computing...
clean all attribute to release memory.
```

It will output two files: `calcuout_RgxDf.tsv.gz` and `calcuout_RgDf.tsv.gz`

#### Help information

The help information of tichr CLI can be accessed by `tichr calcu --help`:

``` text
options:
  -h, --help            show this help message and exit

Input argument:
  readFileList          A list of input files for epigenomic data such as ChIP-
                        seq, ATAC-seq, or CUT&Tag.\ Supported formats include BAM,
                        BigWig, and BedGraph. Multiple files should be provided as
                        a list, e.g.,
                        ["testdata/DNase_rep1.bam","testdata/DNase_rep2.bam"]
  candidateSite         Candidate regulatory sites, this could be a BED3 file
                        (e.g, "testdata/candidatePeak.bed"), or a predefined type
                        specified by a string:
                        "denovo_peak","surronding_bin","onlypromoter".
  candidateGeneFile     A tab-separated file listing candidate genes. It must
                        contain at least six columns in the following order:
                        [chromosome, start, end, gene symbol, gene ID, strand
                        (+/-)]
  gtfile                A tab-separated genome table file, with column 1
                        specifying chromosome names and column 2 indicating
                        chromosome lengths.
  --refGeneFile REFGENEFILE
                        Reference gene file in the same format as
                        candidateGeneFile. This file is used to define gene
                        promoter regions. Typically, it can be the same as
                        candidateGeneFile
  --TSSrange TSSRANGE   Defines the promoter region as transcription start site
                        (TSS) ± this range.
  --S2Gmax S2GMAX       Maximum distance (in base pairs) allowed between a peak
                        and a gene for linking
  --hicfilepath HICFILEPATH
                        Path to hic files in Juicer .hic format.
  --readFileList2 READFILELIST2
                        A second set of epigenome data files, in the same format
                        as readFileList. For example, DNase signals can be
                        provided in readFileList and H3K27ac signals in
                        readFileList2. These two signals are combined using the
                        geometric mean

Process epigenome data arguments:
  --macs2species MACS2SPECIES
                        Used only when candidateSite is set to "denovo_peak".
                        Specifies the effective genome size for MACS2 peak
                        calling. It can be a numeric value (e.g., 1000000000) or a
                        shortcut string (‘hs’ for human, ‘mm’ for mouse, ‘ce’ for
                        C. elegans, ‘dm’ for Drosophila)
  --binResolution BINRESOLUTION
                        Used only when candidateSite is set to “surronding_bin”.
                        Defines the bin size (in base pairs) for creating windowed
                        candidate sites.
  --blackregion BLACKREGION
                        Regions to exclude, such as ENCODE blacklist sites,
                        provided in Bed3 format
  --fixPeakWidth        : Applicable only for given BED3 candidate sites. If set
                        to True, each peak’s width is fixed to 500 bp by centering
                        and extending ±250 b
  --tmpdir TMPDIR       Temporary directory name for intermediate files. Default
                        is a randomly generated name like tichr_tmp_rsDuchihKJ
  --coverageMethod COVERAGEMETHOD
                        Method used to compute coverage. For most users,
                        “coverageBed” is recommended.
  --spmr                Whether to normalize signal by total mapped reads (Signal
                        Per Million Reads). Set to True if you plan to compare Rg
                        or RgX across samples
  --multiBAMmerge MULTIBAMMERGE
                        Strategy to merge multiple replicates. Options: 'mean'
                        (default) or 'sum'.
  --file_type FILE_TYPE
                        Format of epigenomic input files. Supported types: “bam”,
                        “bigwig”, or “bedGraph”

Process Hi-C data arguments:
  --hicRes HICRES       Resolution for Hi-C contact, for example, 10000 for 10kb
                        resolution
  --hicDataType HICDATATYPE
                        could be rawhic_sparse (recommended), matrix_dense (dense
                        matrix for each chromosome), or rawhic_dense (used for
                        ‘strange’ hic files such as that generated by juicertools
                        >2.0. This is the last choice if there are any bugs for
                        the rawhic_sparse mode)
  --hicNormType HICNORMTYPE
                        Normalization type for Hi-C data. Options: 'KR' (Knight-
                        Ruiz), 'VC' (vanilla coverage), 'VC+S' (vanilla coverage +
                        sparse), 'none' (no normalization).
  --juicertool JUICERTOOL
                        Path to juicer_tools.jar file. Only for
                        hicDataType=rawhic_dense. Give a user-difined juicertools
                        jar file to process the hic files.
  --threads THREADS     Number of threads to use for processing Hi-C data. Default
                        is 1.
  --contactNorm CONTACTNORM
                        default: default normalize; abc: similar normalization to
                        the ABC model; oe: observed/expected normalize; 0to1:
                        divide by 95 quantile values; total: divide by the sum of
                        all values, then muliply 1e7.
  --ifUseHiCRef         If set, uses the Hi-C reference file to calculate RgX.

Calculation arguments:
  --weightType WEIGHTTYPE
                        Determines how the site-to-gene weight is calculated.
                        Options: 'hic' (based on Hi-C contact frequency) or
                        'fixed_function' (based on genomic distance).
  --fixedFunctionType FIXEDFUNCTIONTYPE
                        Specifies the function used if
                        weightType=”fixed_function”. Options include: “Sigmoid”,
                        “Exponential”, “Powerlaw”, “NormPL”, “Linear”, “Constant”,
                        “Closest”, or “OnlyPromoter”.
  --halfDistance HALFDISTANCE
                        Distance (in bp) at which the weight decays to 0.5 for
                        supported functions [sigmoid,exponential,powerlaw,linear-
                        half]
  --setpromoter1        If set, sets the RgX ratio of promoter regions to 1
  --threadscalcu THREADSCALCU
                        Not recommended. Number of threads for calculation

Output argument:
  --outname OUTNAME     Output name prefix
  --outhead             Define if the output file with column name
  --noise_ratio NOISE_RATIO
                        If the proportion of the Rgx value to the Rg value of the
                        gene is less than this value, the Rgx value will be set to
                        0.
  --noise_quantile NOISE_QUANTILE
                        If the percentile of the Rgx value among all Rgx values of
                        the gene is less than this value, the Rgx value will be
                        set to 0.
```

#### Full parameters

##### Input argument

-  `readFileList`: A list of input files for epigenomic data such as ChIP-seq, ATAC-seq, or CUT&Tag. Supported formats include BAM, BigWig, and BedGraph. Multiple files should be provided as a list, e.g., ["testdata/DNase_rep1.bam","testdata/DNase_rep2.bam"]
-  `candidateSite`: Candidate regulatory sites, this could be a BED3 file (e.g, "testdata/candidatePeak.bed"), or a predefined type specified by a string: "denovo_peak","surronding_bin","onlypromoter".
-  `candidateGeneFile`: A tab-separated file listing candidate genes. It must contain at least six columns in the following order: [chromosome, start, end, gene symbol, gene ID, strand (+/-)]
-  `gtfile`: A tab-separated genome table file, with column 1 specifying chromosome names and column 2 indicating chromosome lengths.

-   `--refGeneFile`: Reference gene file in the same format as `candidateGeneFile`. This file is used to define gene promoter regions. Typically, it can be the same as `candidateGeneFile`
-   `--TSSrange` : Defines the promoter region as transcription start site (TSS) ± this range.
-   `--S2Gmax`: Maximum distance (in base pairs) allowed between a peak and a gene for linking.
-   `--hicfilepath`: Path to hic files in Juicer .hic format. It is necessary when you are calculating weights by .hic file.
-   `--readFileList2`: A second set of epigenome data files, in the same format as `readFileList`. For example, DNase signals can be provided in `readFileList` and H3K27ac signals in readFileList2. These two signals are combined using the geometric mean


##### Process epigenome data arguments

- `--macs2species`: Used only when `candidateSite`is set to `"denovo_peak"`. Specifies the effective genome size for MACS2 peak calling. It can be a numeric value (e.g., `1000000000`) or a shortcut string (`"hs"`for human, `"mm"`for mouse, `"ce"`for C. elegans, `"dm"`for Drosophila).
- `--binResolution`: Used only when `candidateSite`is set to `"surronding_bin"`. Defines the bin size (in base pairs) for creating windowed candidate sites.
- `--blackregion`: Regions to exclude, such as ENCODE blacklist sites, provided in Bed3 format.
- `--fixPeakWidth`: Applicable only for given BED3 candidate sites. If set to `True`, each peak's width is fixed to 500 bp by centering and extending ±250 bp.
- `--tmpdir`: Temporary directory name for intermediate files. Default is a randomly generated name like `tichr_tmp_rsDuchihKJ`.
- `--coverageMethod`: Method used to compute coverage. For most users, `"coverageBed"`is recommended.
- `--spmr`: Whether to normalize signal by total mapped reads (Signal Per Million Reads). Set to `True`if you plan to compare Rg or RgX across samples.
- `--multiBAMmerge`: Strategy to merge multiple replicates. Options: `'mean'`(default) or `'sum'`.
- `--file_type`: Format of epigenomic input files. Supported types: `"bam"`, `"bigwig"`, or `"bedGraph"`.
- 

##### Process Hi-C data arguments

**Options:**

- `--hicRes`: Resolution for Hi-C contact, for example, `10000`for 10kb resolution.
- `--hicDataType`: Format of Hi-C data. Options: `"rawhic_sparse"`(recommended), `"matrix_dense"`(dense matrix for each chromosome), or `"rawhic_dense"`(used for 'strange' hic files such as that generated by juicertools >2.0).
- `--hicNormType`: Normalization type for Hi-C data. Options: `'KR'`(Knight-Ruiz), `'VC'`(vanilla coverage),  `'none'`(no normalization), etc
- `--juicertool`: Path to `juicer_tools.jar`file. Only for `hicDataType=rawhic_dense`. Give a user-defined juicertools jar file to process the hic files.
- `--threads`: Number of threads to use for processing Hi-C data. Default is `1`.
- `--contactNorm`: Additional normalization methods. Options: `"default"`(default normalization), `"abc"`(similar to ABC model), `"oe"`(observed/expected), `"0to1"`(divide by 95 quantile values), `"total"`(divide by the sum of all values, then multiply by 1e7).
- `--ifUseHiCRef`: If set, uses the Hi-C reference file to calculate `RgX`.
- 

##### Calculation arguments

- `--weightType`: Determines how the site-to-gene weight is calculated. Options: `'hic'`(based on Hi-C contact frequency) or `'fixed_function'`(based on genomic distance).
- `--fixedFunctionType`: Specifies the function used if `weightType="fixed_function"`. Options include: `"Sigmoid"`, `"Exponential"`, `"Powerlaw"`, `"NormPL"`, `"Linear"`, `"Constant"`, `"Closest"`, or `"OnlyPromoter"`.
- `--halfDistance`: Distance (in bp) at which the weight decays to 0.5 for supported functions [sigmoid, exponential, powerlaw, linear-half].
- `--setpromoter1`: If set, sets the RgX ratio of promoter regions to 1.
- `--threadscalcu`: Number of threads for calculation (not recommended for general use).

##### Output argument

- `--outname`: Output name prefix
- `--outhead`: Define if the output file with column name
- `--noise_ratio`: If the proportion of the Rgx value to the Rg value of the gene is less than this value, the Rgx value will be set to 0.
- `noise_quantile`: If the percentile of the Rgx value among all Rgx values of the gene is less than this value, the Rgx value will be set to 0.




## 2. Adjust S2G regulation

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 10px; color: #000000ff;; background-color: #b5cabaff; border-color: #faebcc;">
<p><strong>Important:</strong></p>
<p>
The adjustment is recommended for S2G-level studies, but NOT recommended for gene-level Rg analyses.
</p>
</div>

#### Description

You can check the full usage by type `tichr adjust -h`

``` text
options:
  -h, --help            show this help message and exit

Adjustment Arguments:
  inputRgx              The input RgX DF file, should be standard output from
                        tichr calcu
  inputRg               The input Rg DF file, should be standard output from tichr
                        calcu
  outdir                The output directory where the adjusted results will be
                        saved.
  --RgXhead             Set it if the RgX files contain header lines.
  --tpmFile TPMFILE     The TPM file contains gene expression values in a single-
                        column format, where each row corresponds to the TPM value
                        of a gene.
  --tpmCols TPMCOLS     The column number in the TPM file that contains the TPM,
                        could be multiple columns. like 1,2
  --tpmHead             If set, the first row of the TPM file is ignored. This is
                        useful when the first row contains column headers.
  --tpmGeneID TPMGENEID
                        The column name in the TPM file that contains gene IDs.
  --rankType RANKTYPE   ranktype could be sumrank or diffrak
  --strucTypeList STRUCTYPELIST
                        A comma-separated list of structure types to be used for
                        adjustment. must be supplied in this way
                        'boundary','tad','loop','stripe','compartmentSame'
  --strucFileList STRUCFILELIST
                        A comma-separated list of files containing structure
                        information. must be supplied in this way 'boundary.bed','
                        tad.bed','loop.bed','stripe.bed','compartmentSame.bed'
  --strucWeightList STRUCWEIGHTLIST
                        A comma-separated list of weights corresponding to each
                        structure type. must be supplied in this way 0.5,1.2,5,2,2
```

#### Typical script

The full version of adjust regulation can be found at [here](1.Compute.md#adjust-for-rgx-and-rg)

``` shell
testdata=./Data/Calculate

rgxfile_raw="test_RgxDf.tsv"
rgfile_raw="test_RgDf.tsv"
outdir="adjusted/"
tpmFile="./Data/Calculate/K562.genes.TPM.txt"

boundary="./Data/Calculate/ENCFF621AIY/tad.boundary.point"
tad="./Data/Calculate/ENCFF621AIY/tad.region"
loop="./Data/Calculate/ENCFF621AIY/loop.bedpe"
stripe="./Data/Calculate/ENCFF621AIY/stripe.bedpe"
compartmentSame="./Data/Calculate/ENCFF621AIY/pc1.compart"


tichr adjust $rgxfile_raw $rgfile_raw $outdir \
        --tpmFile $tpmFile --tpmCols 3,4 --tpmGeneID 2 --tpmHead \
        --rankType "sumrank" \
        --strucTypeList "boundary","tad","loop","stripe","compartmentSame"\
        --strucFileList $boundary,$tad,$loop,$stripe,$compartmentSame \
        --strucWeightList 0.5,1.2,5,2,2
```

The output information:

``` text
Start adjust RgX and Rg...
strucTypeList: ['boundary', 'tad', 'loop', 'stripe', 'compartmentSame']
strucFileList: ['./Data/Calculate/ENCFF621AIY/tad.boundary.point', './Data/Calculate/ENCFF621AIY/tad.region', './Data/Calculate/ENCFF621AIY/loop.bedpe', './Data/Calculate/ENCFF621AIY/stripe.bedpe', './Data/Calculate/ENCFF621AIY/pc1.compart']
strucWeightList: [0.5, 1.2, 5.0, 2.0, 2.0]
***Performing adjRaw
***Performing adjStruc
......strucFile, strucType and strucWeight are all list
***Performing adjrank
......Finish the adjustment for Rg and RgX
***Saving adjusted output
Finish adjust RgX and Rg...
```

If you only want to adjust TPM
``` shell
tichr adjust $rgxfile_raw $rgfile_raw $outdir \
        --tpmFile $tpmFile --tpmCols 3,4 --tpmGeneID 2 --tpmHead \
        --rankType "sumrank"
```

#### Full parameter

##### Input and output
- `inputRgx`: The input RgX DF file, should be standard output from tichr calcu.
- `inputRg`: The input Rg DF file, should be standard output from tichr calcu.
- `outdir`: The output directory where the adjusted results will be saved.
- `--RgXhead`: Set it if the RgX files contain header lines.

##### Optioanl 
- `--tpmFile`: The TPM file contains gene expression values in a single-column format, where each row corresponds to the TPM value of a gene.
- `--tpmCols`: The column number(s) in the TPM file that contains the TPM values. Can be multiple columns, e.g., `1,2`.
- `--tpmHead`: If set, the first row of the TPM file is ignored. This is useful when the first row contains column headers.
- `--tpmGeneID`: The column name in the TPM file that contains gene IDs.
- `--ranktype`: Ranking type for adjustment. Options: `'sumrank'`or `'diffrak'`.

- `--strucTypeList`: A comma-separated list of structure types to be used for adjustment. Must be supplied as: `'boundary','tad','loop','stripe','compartmentSame'`.
- `--strucFileList`: A comma-separated list of files containing structure information. Must be supplied as: `'boundary.bed','tad.bed','loop.bed','stripe.bed','compartmentSame.bed'`.
- `--strucWeightList`: A comma-separated list of weights corresponding to each structure type. Must be supplied as: `0.5,1.2,5,2,2`.




## 3.Identify negative regulations

### Description
> Please check more information from [here](5.Negative.md)

The usage of **Tichr**  `negative` command-line tool includes the following parameters:

`tichr neg -h`

```python
usage: __main__.py neg [-h] [--corrtype CORRTYPE] [--minRgx MINRGX]
                       [--minRgxQuantile MINRGXQUANTILE] [--minRgxRatio MINRGXRATIO]
                       [--geneFC_cutoff GENEFC_CUTOFF] [--geneFDR_cutoff GENEFDR_CUTOFF]
                       [--rgxFC_cutoff RGXFC_CUTOFF]
                       rgCtrl rgTreat rgxCtrl rgxTreat outdir outname

options:
  -h, --help            show this help message and exit

Necessary Arguments:
  rgCtrl                Direction to control group Rg file.
  rgTreat               Direction to treated group Rg file.
  rgxCtrl               Direction to control group Rgx file.
  rgxTreat              Direction to treated group Rgx file.
  outdir                Direction for output files.
  outname               Output file prefix.

Options:
  --corrtype CORRTYPE   Could be 'spearman' or 'pearson' (Default: pearson)
  --minRgx MINRGX       Filter the site-to-gene links by RgX value > min Rgx. (Default: 0)
  --minRgxQuantile MINRGXQUANTILE
                        If the percentile of the Rgx value among all Rgx values of the gene is less
                        than this value, the Rgx value will be set to 0.
  --minRgxRatio MINRGXRATIO
                        Filter the site-to-gene links by RgX Ratio > min RgxRatio (Default: 0)
  --geneFC_cutoff GENEFC_CUTOFF
                        Cutoff for the absolute log fold-change (|logFC|) of gene expression.
                        Default=0.5
  --geneFDR_cutoff GENEFDR_CUTOFF
                        Cutoff for the FDR of gene expression. Default=0.05
  --rgxFC_cutoff RGXFC_CUTOFF
                        Cutoff for the absolute log fold-change (|logFC|) of gene-level regulation.
```

### Example usage

Here is an example to show you how to apply Tichr `neg` CLI in a example dataset:

```shell
export MPLBACKEND=Agg
datadir="./Data/Negative"

rgCtrl="${datadir}/Ctrl_RgDf.tsv"
rgTreat="${datadir}/Treat_RgDf.tsv"
rgxCtrl="${datadir}/Ctrl_RgxDf.tsv"
rgxTreat="${datadir}/Treat_RgxDf.tsv"
outdir="NegativeResults"
outname="ERnegative"

python -m tichr neg \
        $rgCtrl $rgTreat $rgxCtrl $rgxTreat \
        $outdir $outname \
        --corrtype pearson --minRgxQuantile 0.01 \
        --minRgxRatio 0.01 --geneFC_cutoff 0.5 \
        --geneFDR_cutoff 0.05 --rgxFC_cutoff 0.5
```

The `neg` method of **Tichr** CLI is consistent with that of the **Tichr** API `negative`. If you are interested in detailed information about certain parameters, you can find more help <a href="https://tichr.readthedocs.io/en/latest/5.Negative.html">here</a>. 

### Output

The success output information

``` text
Merging data frames...
Computing negative regulation...
----------Iteration 0 --------------
----------Iteration 1 --------------
----------Iteration 2 --------------
pearson correlation - Before: 0.0740, After: 0.2498, Difference: 0.1759
ERnegative
Ratio of improved genes :  1.0
Ratio of negative site-to-genes: 0.07304116865869854
Number of genes with negative sites: 116
Percentage of genes  with negative sites: 0.23529411764705882
Number of negative sites: 372
Percentage of negative sites for degs: 0.08391608391608392
------------END------------
---------------------------
---------------------------
---------------------------
---------------------------
```

