# Pipeline with pFind

Download pFind from http://pfind.ict.ac.cn/software/pFind/



## 1. Identification





## 2. Preparation of input file

Extract spectra info from pFind searching results and 

```bash
mkdir pfd
cd pfd
python ../filter_pFind_results.py ../results/IAtev_light.txt ../results/pFind-Filtered_light.spectra > pFind_light.spectra
python ../filter_pFind_results.py ../results/IAtev_heavy.txt ../results/pFind-Filtered_heavy.spectra > pFind_heavy.spectra
python /path/of/CIMAGE/python/pFind2CimageTab.py 20181026_1TO1 "IAtev-L[C]:521.30755:*" "IAtev-H[C]:527.321364:*" "Oxidation[M]:15.9949:#|Carbamidomethyl[C]:57.02146:-|Acetyl[ProteinN-term]:42.010565:@"
```



## 3. Quantification with CIMAGE



```
R --vanilla --args /path/of/IAtev/cimage.params.IAtev 20181026_1TO1 < /path/of/CIMAGE/R/findMs1AcrossSetsFromDTASelect.R 
cd ..
cimage_combine pfd
```

