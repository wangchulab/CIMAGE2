# Pipeline with Mascot



## Identification


## Preparation of input file


## Quantification with CIMAGE
SILAC: cimage SILAC your-cimage-params-file PREFIX
ABPP: cimage_mascot 1mod paramfile modify_AA light_number heavy_number PREFIX 
      cimage_mascot 2mod parameter_file modify_AA1 light_number1 heavy_number1 modify_AA2 light_number2 heavy_number2 PREFIX
PS: modify_AA represents the labelled amino acid (e.g., the modify_AA should be “C” when labelling cysteine). light_number and heavy_number are the number allocated to light and heavy modifications. If you probe can label two different amino acids, use the “2mod” command.

