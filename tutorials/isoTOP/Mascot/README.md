# Pipeline with Mascot

**Note:** The pipeline was only tested on a very old version 2.3.2



## Data preparation

MSConvert is required for converting raw file to mgf format (peak picking with msLevel>1 and remove zero samples can reduce the size of output mgf file significantly).



## Identification

Search with two variable modifications on CYS (light: 521.3074 and heavy: 527.3212) and fixed Carbamidomethyl modification (57.02146).

Put the output **20181026_1TO1_1.csv** into **dta** folder, which is at the same level of mzXML file.




## Quantification with CIMAGE
Open the identification csv file, check the number of modifications, for example

```bash
"Variable modifications","--------------------------------------------------------"

"Identifier","Name","Delta","Neutral loss(es)"
1,"Carbamidomethyl (C)",57.021464
2,"Oxidation (M)",15.994915,0,63.998285
3,"TEV-IA-heavy (C)",527.321225
4,"TEV-IA-light (C)",521.307416
```

Then use **cimage_mascot** for quantification:

```
cimage_mascot 1mod /path/of/cimage.params C 4 3 20181026_1TO1
```

where "**C**" means modified amino acid, "**4**" is number of light probe, "**3**" is number of heavy probe.

