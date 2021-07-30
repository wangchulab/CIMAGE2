# Pipeline with MSFragger

Download MaxQuant from https://www.maxquant.org/



## Identification

Search with light and heavy labels as variable modification separately and search twice, rename the evidence.txt file and put them into **dta** directory

```bash
$>tree dta
├── evidence_light.txt
└── evidence_heavy.txt
```

 

## Quantification

Run **cimage_andromeda** to generate all quantification results

```bash
cd dta
cimage_andromeda /path/of/cimage.params.SILAC no_mod 20181207_SILAC_1LH
cd ..
cimage_combine dta
```

