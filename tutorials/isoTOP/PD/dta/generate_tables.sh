python ../../../../python/csv2CimageTab.py -t HL_ISO_123 \
	-f ../Homo_sapiens_SwissProt_TaxID_9606.fasta \
	-l "IAAL (12C6) com. azo:333:*" \
	-h "IAAH (13C6) com. azo:339:*" \
	-n "M:16:#|Met-loss:0:%|Met-loss+Acetyl:0:%|Acetyl:0:%" \
	-L ../light/HEK_Azo_LtoH_1to1_Test_Whole_Proteome_PSM_Light.txt \
	-H ../heavy/HEK_Azo_LtoH_1to1_Test_Whole_Proteome_PSM_Heavy.txt \
	-j 1
