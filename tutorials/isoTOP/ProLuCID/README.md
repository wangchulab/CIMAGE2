# Pipeline with ProLuCID

Download [ProLuCID](http://fields.scripps.edu/yates/wp/?page_id=821)

Convert raw file to **ms2** and **mzXML** with [Rawconverter](http://fields.scripps.edu/rawconv/)

<img src="../../pics/image-20210821200922153.png" style="zoom:50%;" />



## 1. Identification

Copy or link ms2 file into light and heavy directory, search with light or heavy params

<img src="../../pics/image-20210821202249175.png" alt="image-20210821202249175" style="zoom:80%;" />

```bash
java -Xmx4G -jar /path/of/ProLuCID1_3.jar 20181026_1TO1_1.ms2 search.xml 4
```

Run DTASelect for FDR control of the resulting sqt file

```bash
DTASelect --trpstat --modstat -p 1 -y 2 -m 0 -l Keratin --fp 0.01
```

<img src="../../pics/image-20210821202708580.png" alt="image-20210821202708580" style="zoom:80%;" />



## 2. Preparation of input file

Copy the two DTASelect-filter.txt files into "dta" folder which is located at the same level of 20181026_1TO1_1.mzXML

```bash
mkdir dta
cp /path/of/light/DTASelect-filter.txt dta/DTASelect-filter_20181026_1TO1_light.txt
cp /path/of/heavy/DTASelect-filter.txt dta/DTASelect-filter_20181026_1TO1_heavy.txt
```

<img src="../../pics/image-20210821203026740.png" alt="image-20210821203026740" style="zoom:80%;" />

## 3. Quantification with CIMAGE

Then execute CIMAGE in the "dta" folder, and generate html for visualization

```bash
cd dta
cimage /path/of/cimage.params.IAtev 20181026_1TO1
```

<img src="../../pics/image-20210821203648339.png" alt="image-20210821203648339" style="zoom:80%;" />

```bash
cd ..
cimage_combine dta
```

<img src="../../pics/image-20210821203903883.png" alt="image-20210821203903883" style="zoom:80%;" />



