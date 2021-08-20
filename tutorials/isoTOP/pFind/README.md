# Pipeline with pFind

Download pFind from http://pfind.ict.ac.cn/software/pFind/



## Identification

1. Add the light and heavy tags as two new modifications "IAtev-L[C]" and "IAtev-H[C]"

2. Search with light or heavy tag as variable modification separately

3. Copy **pFind-Filtered.spectra** and extract **all_results.txt** of spectra into **results** directory, rename files accordingly

```
>tree results
results
├── IAtev_heavy.txt
├── IAtev_light.txt
├── pFind-Filtered_heavy.spectra
└── pFind-Filtered_light.spectra
```



## Preparation of input file

Extract spectra info from pFind searching results and generate CIMAGE input tables

```bash
mkdir dta
cd dta
python /path/of/CIMAGE/python/filter_pFind_results.py ../results/IAtev_light.txt ../results/pFind-Filtered_light.spectra > pFind_light.spectra
python /path/of/CIMAGE/python/filter_pFind_results.py ../results/IAtev_heavy.txt ../results/pFind-Filtered_heavy.spectra > pFind_heavy.spectra
python /path/of/CIMAGE/python/pFind2CimageTab.py 20181026_1TO1 "IAtev-L[C]:521.30755:*" "IAtev-H[C]:527.321364:*" "Oxidation[M]:15.9949:#|Carbamidomethyl[C]:57.02146:-|Acetyl[ProteinN-term]:42.010565:@"
```



## Quantification with CIMAGE

Run CIMAGE and generate output webpage

```
R --vanilla --args /path/of/IAtev/cimage.params.IAtev 20181026_1TO1 < /path/of/CIMAGE/R/findMs1AcrossSetsFromDTASelect.R 
cd ..
cimage_combine dta
```

