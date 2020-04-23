# CIMAGE_Lite
Lite version of CIMAGE, scripts for quantification only. (no web server)
Contact: chuwang@pku.edu.cn or wendao@pku.edu.cn

## Requirements

R (tested on 3.4.1 and 3.4.4)

```
sudo apt install r-base r-cran-ncdf4
```

Python (tested on version 2.7)

Perl (tested on perl5 version 18 and 26)

Other libray

```
sudo apt install libxml2-dev libnetcdf-dev zlib1g-dev
```

XCMS (tested on version 1.51.1 and 1.52.0)

```
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("xcms"))
```

## Installation

Add the path of CIMAGE to environment variable: **CIMAGE_PATH**

Add the path of CIMAGE/Shell to environment variable: **PATH**

## Useage

### Input:
A. For LC-MS/MS data, CIMAGE takes .mzXML files as input, which can be easily converted from raw files by variouty of tools like [ReAdW](http://tools.proteomecenter.org/wiki/index.php?title=Software:ReAdW) and [RawConverter](http://fields.scripps.edu/rawconv/).

B. Before CIMAGE quantification, database search should be performed first and the output file is then used as the input for CIMAGE. CIMAGE can be connected with variety of database search engines (ProLuCID, Mascot, Andromeda, pFind and MSFragger).

### Analysis workflow:

The demo data and results can be downloaded from the following link: https://drive.google.com/drive/folders/1yLNTHNciC9vllmzdA7vuPHmrxjjQ04ce?usp=sharing

A. Make a folder such as "isoTOP-ABPP", upload LC-MS/MS data in mzXML format and create a folder (e.g. dta).

B. Before quantification, you need to edit your cimage.params file to point to the right light/heavy chemical composition files. Template cimage.params file and common light/heavy chemical composition files can be found in the params folder:

```
IA _tev – isoTOP-ABPP with IA-alkyne probe & Tev tag
SILAC – SILAC with Arg10 & Lys8
```

Light/heavy tables are strictly tab-delimited text files, so it is better to use EXCEL or notepad++ to import, edit and then export it.

C. Enter the “dta” folder and upload the database search results. In database search process, the LC-MS/MS data should be searched for light and heavy labelled peptides separately to generate two results containing either light or heavy peptides. By default, ProLuCID is recommended as it is flexible to fit different types of quantitative strategies, but you can also feel free to use your comfortable one.
