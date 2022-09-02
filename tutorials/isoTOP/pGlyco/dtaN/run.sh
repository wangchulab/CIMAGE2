python ../../../../python/pGlyco2CimageTab.py "LPG:0:*" "HPG:0:*" "Oxidation[M]:0:#|Carbamidomethyl[C]:0:-" ../N-light ../N-heavy > log.txt

for d in $(grep "G:" log.txt | cut -d" " -f2)
do
cd $d
python ../../../../../python/generate_cimage_table.py tab.txt
cd ..
done

grep "G:" log.txt | awk '{print "cd", $2"; R --vanilla --args ../../cimage.params.PG.txt Lumos3_HomeMadeCol_ChenXing_LJL_Hela_PGHL_75um_50cm_360min_HCD30_EThcDSA35_20220106_F1 < ~/share/CIMAGE_Lite/R/findMs1AcrossSetsFromDTASelect_v2.R; cd .." }' > jobs.list
