# Pipeline with MSFragger

Download MSFragger from https://github.com/Nesvilab/MSFragger

Download philosopher from https://philosopher.nesvilab.org/



## Identification

1. Convert raw file to mzML with [MSConvert](https://sourceforge.net/projects/proteowizard/), create a **dta** directory, link mzML and copy fasta (with reversed decoys) into it (there will be two copies of mzML file, one for identification and one for quantification)

2. Search with two (light and heavy) variable modifications, see **closed_fragger.params**, you need to delete symbol "**'**"
from the description line within fasta file, it conflicts with html output format
   ```bash
   $>tree dta
   ├── 20181026_1TO1_1.mzML -> ../20181026_1TO1_1.mzML
   ├── closed_fragger.params
   └── target.fasta -> /path/of/fasta
   ```

3. Run MSFragger and philosopher

```
cd dta
java -Xmx32G -jar /path/of/MSFragger-2.4.jar closed_fragger.params 20181026_1TO1_1.mzML
philosopher workspace --init
philosopher database --annotate target.fasta --prefix Reverse_
philosopher peptideprophet --nonparam --expectscore --decoyprobs --ppm --accmass --decoy Reverse_ --database target.fasta *.pepXML
philosopher proteinprophet --maxppmdiff 2000000 --output combined *.pep.xml
philosopher filter --sequential --razor --mapmods --prot 0.01 --tag Reverse_ --pepxml . --protxml combined.prot.xml
philosopher report --mzid
philosopher workspace --clean
```

4. **psm.tsv** file is expected for the next step



## Preparation of input file

Generate input tables for CIMAGE

```bash
/path/of/CIMAGE/python/pepPSMs2CimageTab.py 20181026_1TO1 target.fasta 'C:464.28595:*' 'C:470.29976:*' 'M:15.99490:#|n:42.01060:@' . .
```
The options are

- name of raw file

- fasta file

- modifications of light probe

- modifications of heavy probe

- modifications of normal peptide

- location of psm.tsv for light label

- location of psm.tsv for heavy label (same as light if you only search once)

  

## Quantification with CIMAGE

Run CIMAGE (it supports both mzXML or mzML)

```
R --vanilla --args /path/of/IAtev/cimage.params.IAtev 20181026_1TO1 < /path/of/CIMAGE/R/findMs1AcrossSetsFromDTASelect_v2.R #for mzML
cd ..
cimage_combine dta
```

