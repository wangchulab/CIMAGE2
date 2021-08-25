# Table Generation

If you want to add your own modifications for CIMAGE, first you need to know the **elemental composition** for each modification, then write a txt file specifying the mark you used in **search.xml** (or search engines other than ProLuCID).

    ISO-L-323 #:O,1 @:C,2:H,2:O,1 *:C,16:H,29:N,5:O,2
    ISO-H-340 #:O,1 @:C,2:H,2:O,1 *:C,16:H,30:N,5:O,3

Each line generates one table, the first column is the file name of the table, and the following column define the elemental composition of each modification and mark, then type

```
python /path/of/CIMAGE/python/generate_cimage_table.py cy2.txt
```

Two tables will be generated and can be validated by

```
$> python /path/of/CIMAGE/python/check_cimage_table.py ISO-L-323.table.txt
* 323.232125
# 15.994915
@ 42.010565
$> python /path/of/CIMAGE/python/check_cimage_table.py ISO-H-340.table.txt
* 340.234865
# 15.994915
@ 42.010565
```



Notes: tables should match with search params exactly.

