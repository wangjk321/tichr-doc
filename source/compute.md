# Compute Rg and RgX

This page provide necessary steps for obtain Rg for each gene, and Rgx for each site-to-gene links.


## A quick start

``` python
RPobj = Tichr(candidateSite,[bam,],gtfile,genebedfile,refgene_file=genebedfile,peakToGeneMaxDistance=5000000)
RPobj.makeSiteBed()
RPobj.makeSiteBdg("coverageBed",spmr = True)
RPobj.hicfilepath = hic_file
RPobj.proceessHiC(25000,'rawhic_sparse',"VC_SQRT",
                    threads=12,further_normalize_type='no_further')
RPobj.computeAllGene("hic")
RPobj.RgxDf.to_csv("resultdf/"+name+"_RP_RgxDf.tsv",header=None,sep="\t",index=None)
RPobj.RgDF.to_csv("resultdf/"+name+"_RP_RgDf.tsv",header=None,sep="\t",index=None)

```

## Make Tichr object

``` python
candidatesite="testdata/candidatePeak.tsv"
candidateGeneFile="testdata/candidateGene.tsv"
DNase_rep1="testdata/DNase_rep1.bam"
DNase_rep2="testdata/DNase_rep2.bam"
readFileList=[DNase_rep1,DNase_rep2,]
gtfile="testdata/hg38_gt.tsv"

refgene_file="testdata/all_coding_gene.tsv"
H3K27ac="testdata/H3K27ac.bam"
readFileList2=[H3K27ac,]
hicfilepath="testdata/ENCFF621AIY.hic"

tichrobj = Tichr(candidatesite,readFileList,gtfile,candidateGeneFile,refgene_file=refgene_file,
                 ifTSSrange=500,peakToGeneMaxDistance=100000,
                 hicfilepath=hicfilepath,readFileList2=readFileList2)

```

**Required parameters**

- `candidatesite`: candidate regulatory sites, this could be `path-to-bed3-file` (e.g, `"testdata/candidatePeak.tsv"`), or one of the fixed sites by assign the strings `"denovo_peak"`,`"surronding_bin"`,`"onlypromoter"`. type=str,default="onlypromoter"
- `readFileList`: list of files for epigenome data from ChIP-seq,ATAC-seq,CUT&TAG etc.,the format could be bam / bigwig /bedGraph. multiple files should be supplied by a list such as `["testdata/DNase_rep1.bam","testdata/DNase_rep2.bam"]`
- `gtfile`: Tab-separated genome table file, column1: chromosome name, column2: chromosome length.
- `candidateGeneFile`: candidate gene files, tab-sperated, at least six columns: [chr,start,end,gene symbol,gene ID,+or- strands]

**Optional parameters**
- `refgene_file`: all reference gene files, same format as candidateGene. Usually you can specify it the same file for `candidateGeneFile`, it only used to define gene promoters. 
- `ifTSSrange`: Range of promoters; i.e., transcription start sites ± TSSrange,type=int,default=500.
- `peakToGeneMaxDistance`: "Max distance (bp) for peak-to-gene links",type=int,default=100000
- `hicfilepath`: path for hic files in Juicer format or matrix_dense
- `readFileList2`: Second types of epigenome data, should be the same format as epigenomeData 1. For example, you can use DNase signal in readFileList and H3K27ac signal in readFileList2. Two types of signal are combined by geometric mean

## Generate candidate sites

``` python
encodeblack="testdata/hg38-blacklist.v2.bed"
tichrobj.makeSiteBed(macs2species='hs',binResolution=100,
                blackregion=encodeblack,tmpdir=None,fixPeakWidth=False)
```

- `macs2species`: Only for `denovo_peak` candidateSite. Effective genome size of macs2. It can be number (1000000000), or shortcuts ('hs','mm','ce', 'dm'), Default:hs''',type=str,default="hs"
- `binResolution`: Only for `surronding_bin`  candidateSite. Bin sizes(bp) for make windowed sites. default: 100
- `fixPeakWidth`: Only for given candidate sites. If True, fix the peak width by center +- 250 bp, generating a width=500 peaks.
- `blackregion`: Regions that should be removed such as ENCODE black list. Bed3 format.
- `tmpdir`: temp directory name. default is by random characters such as tichr_tmp_rsDuchihKJ

You can check the obtained sites file by the attibute `tichrobj.candidatesite_file`

## Calculate epigenome activities

``` python
tichrobj.makeSiteBdg("coverageBed",spmr = False,multiBAMmerge='mean',file_type="bam")
```

- `coverageMethod`: For most users, use "coverageBed"
- `spmr`: If normalize signal by the total reads. If you are comparing the Rg or RgX among samples, it is recommend to set as True.
- `multiBAMmerge`: Could be 'mean' or 'sum'. default=mean
- `file_type`: file type for the epigenomes, could be [bam, bigwig, bedGraph]

Addtional parameter (usually not used)
- quantileref=None
- quantile_method=None
- signaltype=None
- separatepromoter=False
- signal2type=None

The main obtained attributte is `tichrobj.candidatesite_coverage`

| Chromosome | Start Position | End Position | Value       | ifPromoter |
|------------|----------------|--------------|-------------|-------|
| chr1       | 25933091       | 25934590     | 341.499634  | 0     |
| chr1       | 26378064       | 26378863     | 552.494344  | 0     |
| chr1       | 26508509       | 26508806     | 21.633308   | 1     |
| chr1       | 27867080       | 27867532     | 188.403822  | 0     |


## Process Hi-C files (Optional)
>If you plan to use some functional such as exponential, you can skip this step.

**example usage**
``` python
tichrobj.proceessHiC(50000,'rawhic_sparse',"VC_SQRT",
                   threads=12,further_normalize_type='no_further')
```

**all parameters**
``` python
proceessHiC(hicRes,hicDataType,hicNormType,juicertool=None,
                threads=8,further_normalize_type='abc'):
```

- `hicRes`: resolution for hic contact matrix
- `hicDataType`: could be `rawhic_sparse` (recommended), `matrix_dense` (dense matrix for each chromosome), or `rawhic_dense` (used for 'strange' hic files such as that generated by juicertools >2.0. This is the last choice if there are any bugs for the rawhic_sparse mode)
- `hicNormType`: could be [VC_SQRT, KR, SCALE, VC]
- `juicertool`: Only for `hicDataType=rawhic_dense`. Give a user-difined juicertools jar file to process the hic files. 
- `threads`: number of threads for processing hic files. This parallel computing is achieved by computing for each chromosomes separately.
- `further_normalize_type`: could be `no_further`, `abc`, `oe`, `0to1`,`total`. 
    - `no_further`: default normalize
    - `abc`: similar normalization to the ABC model
    - `oe`: observed/expected normalize
    - `0to1`: divide by 95% quantile values. 
    - `total`: divide by the sum of all values, then muliply 1e7.

The main obtained information is `tichrobj.nomhicdf`, which provide the contact information for all chr

``` python
tichrobj.nomhicdf["chr19"]
```

| Bin Pair        | posX   | posY   | Counts       |
|-----------------|--------|--------|--------------|
| 50000to50000    | 50000  | 50000  | 2477.239014  |
| 200000to200000  | 200000 | 200000 | 1096.694336  |
| 50000to250000   | 50000  | 250000 | 52.562267    |


## Calculate Rg and Rgx

example usage:
``` python
tichrobj.computeAllGene("hic")
```

full usage:
``` python
tichrobj.computeAllGene(weightType,fixedFunctionType='rp-classic',halfDistance=10000,
                 given_gamma=1.024238616787792, given_scale = 5.9594510043736655,
                 ref_gamma = 0.87, ref_scale = -4.80 + 11.63 * 0.87, hicmindistance=5000,
                 setpromoter1=False,threads=1,ifUseHiCRef=False):
```

- `weightType`: If use a fixed function or Hi-C information to calculate the weight for RgX,could be [hic,fixed_function]
- `fixedFunctionType`: If you choose a fixed function, decide which type. Could be [closest,sigmoid,exponential,powerlaw,normPowerLaw,constant1,linear-half,linear]''',default="exponential"
- `halfDistance`: Parameter for [sigmoid,exponential,powerlaw,linear-half]
- `setpromoter1`: If set the Rgx ratio of Rgx to 1
- `threads`: Not recommended.number of threads for calculation


Addtional parameter (usually not used)
- given_gamma
- given_scale
- ref_gamma
- ref_scale
- hicmindistance
- logRgX=False
- ifUseHiCRef

The main obtained information is `tichrobj.RgxDf`:

| Peak Chr | Peak Start | Peak End | Activity     | Gene Symbol         | Gene Chr | Gene Start | Gene End | Strand | Gene ID             | Weight     | Rgx Raw Value | Rgx Percent |
|----------|------------|----------|--------------|----------------------|----------|------------|----------|--------|----------------------|------------|----------------|-------------|
| chr1     | 25933091   | 25934590 | 341.499634   | ENST00000516882.1    | chr1     | 26006216   | 26006216 | +      | ENSG00000252691.1    | 984.895508 | 336341.5       | 1.0         |
| chr1     | 26378064   | 26378863 | 552.494344   | ENST00000319041.6    | chr1     | 26280122   | 26280122 | +      | ENSG00000142669.9    | 1761.7518  | 973357.9       | 0.077389         |
| chr1     | 26378064   | 26378863 | 552.494344   | ENST00000492808.1    | chr1     | 26317958   | 26317958 | +      | ENSG00000169442.4    | 2596.0479  | 1434302.0      | 0.893021       |

and `tichrobj.RgDF` (the input gene file with the last column represent the Rg values):

| Chromosome | Start     | End       | Symbol       | Gene ID              | Strand | Value         |
|------------|-----------|-----------|----------------------|-----------------------|--------|---------------|
| chr1       | 25272548  | 25272548  | geneA    | ENSG00000187010.14    | +      | 0.000000       |
| chr1       | 26006216  | 26006216  | geneB    | ENSG00000252691.1     | +      | 336341.455414  |
| chr1       | 26183204  | 26183204  | geneC    | ENSG00000142675.13    | +      | 0.000000       |


You can save the main output by 

``` python
tichrobj.RgxDf.to_csv("resultdf_hic/name_RgxDf.tsv.gz",header=None,sep="\t",index=None,compression='gzip')
tichrobj.RgDF.to_csv("resultdf_hic/name_RgDf.tsv.gz",header=None,sep="\t",index=None,compression='gzip')
```

## Manipulate Tichr object

As a object-oriented programing, you can esaily re-assign parameters in the intermediate step and re-run the program.

For example, you want to use a new hic file, you can try:

``` python
tichrobj.hicfilepath = new_hic_file
tichrobj.proceessHiC(50000,'rawhic_sparse',"VC_SQRT")
tichrobj.computeAllGene("hic")
```

Another important thing is when you use super high resolution hic or very large number of peaks-to-genes, the tichr object will consume much memory. To check the memory used, you can try:

``` python
def format_bytes(size):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024

from pympler import asizeof
size_bytes = asizeof.asizeof(tichrobj)
print("Raw size (bytes):", size_bytes)
print("Readable size:", format_bytes(size_bytes))
```
The output is like:
```
Raw size (bytes): 26964769808
Readable size: 25.11 GB
```

To clean the major memory consumption attribute, you can use
``` python
tichrobj.clean()
```
Then you can found the memory is cleaned.
```
Raw size (bytes): 1387440
Readable size: 1.32 MB
```