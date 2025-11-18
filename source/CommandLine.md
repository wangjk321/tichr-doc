# Command Line Interface

The core function of **Tichr** -  `compute` and `negative` -  can be used in the form of a command-line tool. 

Unlike the **Tichr** API, the **Tichr** CLI is more convenient to use. It does not require distributed processing and can directly obtain **Rg** and **RgX** results through a single command line. However, in terms of flexibility, it falls short compared to the **Tichr** API. 



## 1.Compute

The usage of **Tichr**  `calculate` command-line tool includes the following parameters:

``` shell
$ tichr calcu [-h] [--refgene_file REFGENE_FILE] [--TSSrange TSSRANGE]
                         [--peakToGeneMaxDistance PEAKTOGENEMAXDISTANCE] [--hicfilepath HICFILEPATH]
                         [--readFileList2 READFILELIST2] [--outdir OUTDIR] [--macs2species MACS2SPECIES]
                         [--binResolution BINRESOLUTION] [--blackregion BLACKREGION] [--fixPeakWidth FIXPEAKWIDTH]
                         [--tmpdir TMPDIR] [--coverageMethod COVERAGEMETHOD] [--spmr]
                         [--multiBAMmerge MULTIBAMMERGE] [--file_type FILE_TYPE] [--hicRes HICRES]
                         [--hicDataType HICDATATYPE] [--hicNormType HICNORMTYPE] [--juicertool JUICERTOOL]
                         [--threads THREADS] [--further_normalize_type FURTHER_NORMALIZE_TYPE] [--ifUseHiCRef]
                         [--weightType WEIGHTTYPE] [--fixedFunctionType FIXEDFUNCTIONTYPE]
                         [--halfDistance HALFDISTANCE] [--setpromoter1] [--threadscalcu THREADSCALCU]
                         [--tpmfile TPMFILE] [--tmpcolrep TMPCOLREP] [--ignorehead] [--tmpgeneID TMPGENEID]
                         [--structureTypeList STRUCTURETYPELIST] [--structureFileList STRUCTUREFILELIST]
                         [--structureWeightList STRUCTUREWEIGHTLIST] [--ranktype RANKTYPE]
                         readFileList candidateSite candidateGeneFile gtfile
```

Here is an example to show you how to apply Tichr `calcu` CLI in a example dataset:

```shell
testdata=../DataSmall/Calculate
candidatesite="$testdata/candidatePeak.bed"
candidateGeneFile="$testdata/candidateGene.tsv"
DNase_rep1="$testdata/DNase_rep1.bigwig"
DNase_rep2="$testdata/DNase_rep2.bigwig"
H3K27ac="$testdata/H3K27ac.bigwig"
gtfile="$testdata/hg38_gt.tsv"
hicfilepath="$testdata/ENCFF621AIY.hic"
refgene_file="$testdata/all_coding_gene.tsv"
encodeblack="$testdata/hg38-blacklist.v2.bed"  

tichr calcu $DNase_rep1,$DNase_rep2 \
	$candidatesite $candidateGeneFile $gtfile \
	--blackregion $encodeblack --file_type bw --readFileList2 $H3K27ac \
        --hicfilepath $hicfilepath --refgene_file $refgene_file \
        --hicRes 50000 --hicNormType VC_SQRT --threads 12 \
        --further_normalize_type no_further --weightType hic \
        --outdir calcuout
```

The `calcu` method of **Tichr** CLI is consistent with that of the **Tichr** API `compute`. If you are interested in detailed information about certain parameters, you can find more help <a href="https://tichr.readthedocs.io/en/latest/1.Compute.html">here</a>. 

Next, I will introduce the functions corresponding to different types of parameters:



### Input and output argument

**Necessary Arguments:** 

-  `--readFileList`: A list of input files for epigenomic data such as ChIP-seq, ATAC-seq, or CUT&Tag. Supported formats include BAM, BigWig, and BedGraph. Multiple files should be provided as a list, e.g., ["testdata/DNase_rep1.bam","testdata/DNase_rep2.bam"]
-  `--candidateSite`: Candidate regulatory sites, this could be a BED3 file (e.g, "testdata/candidatePeak.bed"), or a predefined type specified by a string: "denovo_peak","surronding_bin","onlypromoter".
-  `--candidateGeneFile`: A tab-separated file listing candidate genes. It must contain at least six columns in the following order: [chromosome, start, end, gene symbol, gene ID, strand (+/-)]
-  `--gtfile`: A tab-separated genome table file, with column 1 specifying chromosome names and column 2 indicating chromosome lengths.

**Options**:

-   `--refgene_file`: Reference gene file in the same format as `candidateGeneFile`. This file is used to define gene promoter regions. Typically, it can be the same as `candidateGeneFile`
-   `--TSSrange` : Defines the promoter region as transcription start site (TSS) ± this range.
-   `--peakToGeneMaxDistance`: Maximum distance (in base pairs) allowed between a peak and a gene for linking.
-   `--hicfilepath`: Path to hic files in Juicer .hic format. It is necessary when you are calculating weights by .hic file.
-   `--readFileList2`: A second set of epigenome data files, in the same format as `readFileList`. For example, DNase signals can be provided in `readFileList` and H3K27ac signals in readFileList2. These two signals are combined using the geometric mean
-   `--outdir`: Results output directory



### Process epigenome data arguments

**Options:**

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

### Process Hi-C data arguments

**Options:**

- `--hicRes`: Resolution for Hi-C contact, for example, `10000`for 10kb resolution.
- `--hicDataType`: Format of Hi-C data. Options: `"rawhic_sparse"`(recommended), `"matrix_dense"`(dense matrix for each chromosome), or `"rawhic_dense"`(used for 'strange' hic files such as that generated by juicertools >2.0).
- `--hicNormType`: Normalization type for Hi-C data. Options: `'KR'`(Knight-Ruiz), `'VC'`(vanilla coverage), `'VC+S'`(vanilla coverage + sparse), `'none'`(no normalization).
- `--juicertool`: Path to `juicer_tools.jar`file. Only for `hicDataType=rawhic_dense`. Give a user-defined juicertools jar file to process the hic files.
- `--threads`: Number of threads to use for processing Hi-C data. Default is `1`.
- `--further_normalize_type`: Additional normalization methods. Options: `"default"`(default normalization), `"abc"`(similar to ABC model), `"oe"`(observed/expected), `"0to1"`(divide by 95 quantile values), `"total"`(divide by the sum of all values, then multiply by 1e7).
- `--ifUseHiCRef`: If set, uses the Hi-C reference file to calculate `RgX`.
- 

### Calculation arguments

**Options:**

- `--weightType`: Determines how the site-to-gene weight is calculated. Options: `'hic'`(based on Hi-C contact frequency) or `'fixed_function'`(based on genomic distance).
- `--fixedFunctionType`: Specifies the function used if `weightType="fixed_function"`. Options include: `"Sigmoid"`, `"Exponential"`, `"Powerlaw"`, `"NormPL"`, `"Linear"`, `"Constant"`, `"Closest"`, or `"OnlyPromoter"`.
- `--halfDistance`: Distance (in bp) at which the weight decays to 0.5 for supported functions [sigmoid, exponential, powerlaw, linear-half].
- `--setpromoter1`: If set, sets the RgX ratio of promoter regions to 1.
- `--threadscalcu`: Number of threads for calculation (not recommended for general use).
- 

### Adjustment arguments

**Options:**

- `--tpmfile`: The TPM file contains gene expression values in a single-column format, where each row corresponds to the TPM value of a gene.
- `--tmpcolrep`: The column number(s) in the TPM file that contains the TPM values. Can be multiple columns, e.g., `1,2`.
- `--ignorehead`: If set, the first row of the TPM file is ignored. This is useful when the first row contains column headers.
- `--tmpgeneID`: The column name in the TPM file that contains gene IDs.
- `--structureTypeList`: A comma-separated list of structure types to be used for adjustment. Must be supplied as: `'boundary','tad','loop','stripe','compartmentSame'`.
- `--structureFileList`: A comma-separated list of files containing structure information. Must be supplied as: `'boundary.bed','tad.bed','loop.bed','stripe.bed','compartmentSame.bed'`.
- `--structureWeightList`: A comma-separated list of weights corresponding to each structure type. Must be supplied as: `0.5,1.2,5,2,2`.
- `--ranktype`: Ranking type for adjustment. Options: `'sumrank'`or `'diffrak'`.



## 2.Negative

The usage of **Tichr**  `negative` command-line tool includes the following parameters:

```python
usage: __main__.py neg [-h] [--outdir OUTDIR] [--rg_ctrl RG_CTRL] [--rg_treat RG_TREAT] [--rgx_ctrl RGX_CTRL] [--rgx_treat RGX_TREAT] [--chrrtype CHRRTYPE] [--min_rgx MIN_RGX] [--min_rgx_ratio MIN_RGX_RATIO]
```

**Necessary Arguments:**
  `--outdir OUTDIR`:  Direction for output files.
  `--rg_ctrl RG_CTRL`:  Direction to control group Rg file.
  `--rg_treat RG_TREAT`: Direction to treated group Rg file.
  `--rgx_ctrl RGX_CTRL`: Direction to control group Rgx file.
  `--rgx_treat RGX_TREAT`: Direction to treated group Rgx file.

**Options:**
  `--corrtype`:  Could be '`spearman`' or '`pearson`' (Default: pearson)
  `--min_rgx`:  Filter the site-to-gene links by RgX value > min Rgx. (Default: 0.1)
  `--min_rgx_ratio`  Filter the site-to-gene links by RgX Ratio > min RgxRatio (Default: 0.01)

Here is an example to show you how to apply Tichr `neg` CLI in a example dataset:

```shell
datadir="../Data/Negative"

rgCtrl_path="${datadir}/Ctrl_RgDf.tsv"
rgTreat_path="${datadir}/Treat_RgDf.tsv"
rgxCtrl_path="${datadir}/Ctrl_RgxDf.tsv"
rgxTreat_path="${datadir}/Treat_RgxDf.tsv"

tichr neg \
--rg_ctrl $rgCtrl_path --rg_treat $rgTreat_path \
--rgx_ctrl $rgxCtrl_path --rgx_treat $rgxTreat_path \
--outdir neg_cli 
```

The `neg` method of **Tichr** CLI is consistent with that of the **Tichr** API `negative`. If you are interested in detailed information about certain parameters, you can find more help <a href="https://tichr.readthedocs.io/en/latest/5.Negative.html">here</a>. 

