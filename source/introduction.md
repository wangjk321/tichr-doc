# Introduction and Install

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #8a6d3b;; background-color: #fcf8e3; border-color: #faebcc;">
<b>TICHR is underdevelopment, this page is only for developer usage. For general users, please wait the release of our formal version.</b>
</div>


## What is TICHR

<img src="_static/logo.png" style="zoom:30%;" />

Analyzing transcriptional regulation at multiple omics levels plays a guiding role in understanding cellular activities and human diseases. However, due to the lack of direct correlations among different omics, there are many challenges in integrating multi-omics data. Here, we developed a bioinformatics tool called TICHR for studying transcriptional regulation by integrating epigenome data (e.g., ChIP-seq), 3D genome data (e.g.,Hi-C), and transcriptome data (e.g., RNA-seq). TICHR defines the regulatory effect of loci to genes as a function of epigenetic modifications and chromosomal structure information, and then uses transcriptomics to assess the degree of transcriptional activation or inhibition caused by this regulation potential, establishing a model for the correspondence and regulatory strength between loci and genes on a genome-wide scale. 

-----

## installation

TICHR is developed based on `Python 3.11.5`. Therotically, any version newer than 3.11.5 will support TICHR. However, if you face any problem, please try to use conda to create a independent enviroment with python 3.11.5 .

### Requirement

Python packages:
- pandas
- numpy
- matplotlib
- pyBigWig
- scipy
- sklearn
- seaborn
- hicstraw
- collections
- concurrent
- functools
- statsmodels
- pyranges

Optional
- rpy2
- joblib
- multiprocessing

### Installations

You can use PyPi to install the latest version of Tichr

``` shell
pip install tichr
```

Alternatively, you can use the github source code to install

``` shell
git clone XXXXXX
cd XXXX
python setup.py install
```

### Usage

There are two ways to use Tichr

1. By command line, this is the easiest way for most users. After installation, you can try:
``` shell
tichr --help
```

2. By python module. You can import the python pacakges like:
``` python
import sys
sys.path.append('/home/wang/github/Tichr-CLI/tichr')
from tichr import *
```