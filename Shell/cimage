#!/bin/bash

echo command excuted $0 $@

if [ $# -lt 1 ]; then
    echo Usage: $0 cimage.params [no_png] set_1 set_2 ...
    exit -1
fi
param=$1
shift;
if [ ! -f "$param" ];
then
    echo cannot find cimage.param file : $param
    exit -1
fi
no_png=$1
if [ "$no_png" == "no_png" ]; then
    shift;
fi

mzxml=$@

ndta=$(\ls DTASelect-filter_*.txt | wc -l)
if [ "$ndta" -eq 0 ];
then
    echo No DTASelect-fitler files found and please double check your working path. Exiting ...
    exit -1
fi

echo -n > tmp.ipi_name
echo -n > tmp.key_scan
echo -n > tmp.seq_mass
# tagged and fwd
date
echo Parsing DTASelect files
for p11 in $(echo $mzxml)
do
    p1=$(echo $p11 | sed 's/_HL$//g')
    echo $p1
    for p2 in $(ls DTASelect-filter_$p1\_*.txt);
    do
	if grep -q "IPI:IPI" $p2; then
	    IPIDB="1"
	else
	    if grep -q " GN=" $p2; then
		IPIDB="2";
	    else
		IPIDB="0"
	    fi
	fi
	HL=$(echo $p2 | sed 's/\.txt$//g' | awk -F "_" '{print $NF}')
	$CIMAGE_PATH/python/tagDTASelect.py $p2 > $p2.tagged;
	cat $p2.tagged | grep ^cimagepep | grep $p1 | sed 's/^cimagepep\-//g' | sed 's/IPI://g' | sed 's/'sp\|'//g' | sed 's/'tr\|'//g' | awk '{print $1}' > $p2.tmp.ipi
	cat $p2.tagged | grep ^cimagepep | grep $p1 | awk '{print $2}' > $p2.tmp.FileName
	cat $p2.tagged | grep ^cimagepep | grep $p1 | awk -v HL=$HL '{print $3, HL}' > $p2.tmp.xcorr
	cat $p2.tagged | grep ^cimagepep | grep $p1 | awk '{print $NF}' > $p2.tmp.peptide
	cat $p2.tmp.peptide | sed 's/\*//g' | sed 's/\#//g' | sed 's/\@//g' | cut -f2 -d \.  > $p2.tmp.sequence
	for p3 in $(cat $p2.tmp.sequence | sort | uniq )
	do
	    echo -n "$p3 "
	    $CIMAGE_PATH/python/peptideCalcMass.py $p3 mono;
	done > $p2.tmp.uniq.mass
	for p3 in $(cat $p2.tmp.sequence)
	do
	    grep "^$p3 " $p2.tmp.uniq.mass | awk '{print $2}'
	done > $p2.tmp.mass
	cat $p2.tmp.FileName | awk -F "." '{print $1}' | awk -F "_" '{print $NF}' > $p2.tmp.segment
	cat $p2.tmp.FileName | awk -F "." '{print $2}'  > $p2.tmp.scan
	cat $p2.tmp.FileName | awk -F "." '{print $NF}' > $p2.tmp.charge
	cat $p2.tmp.FileName | awk -v run=$p1 '{print run}' > $p2.tmp.run
	paste -d":" $p2.tmp.ipi $p2.tmp.peptide $p2.tmp.charge $p2.tmp.segment > $p2.tmp.key
	paste -d " " $p2.tmp.run $p2.tmp.scan $p2.tmp.mass $p2.tmp.xcorr $p2.tmp.key >> tmp.key_scan
	rm -rf $p2.tmp.*
	cat $p2.tagged | grep ^cimageipi | sed 's/^cimageipi\-//g' | sed 's/IPI://g' | sed 's/'sp\|'//g' | sed 's/'tr\|'//g' | awk '{print $1} '| awk -F"|" '{print $1}'> tmp.ipi
	## name deliminator "Gene_Symbol=" or "Full="
	if [ $IPIDB == "1" ]; ## ipi database
	then
	    cat $p2.tagged | grep ^cimageipi | awk -F "\t" '{print $NF}' | awk -F"l=" '{print $NF}'| cut -c1-50 | sed -e s/^\-/_/g > tmp.name
	else
	    if [ $IPIDB == "2" ]; then # standard uniprot database
##	    cat $p2.tagged | grep ^cimageipi | awk -F "\t" '{print $NF}' | awk -F"GN=" '{print $NF}'| awk '{print $1}' > tmp.name_gene
		cat $p2.tagged | grep ^cimageipi | awk -F"GN=" '{print $NF}'| awk '{print $1}' > tmp.name_gene
##	    cat $p2.tagged | grep ^cimageipi | awk -F "\t" '{print $NF}' | awk -F"OS=" '{print $1}'| cut -c1-50 | sed -e s/^\-/_/g > tmp.name_desc
		cat $p2.tagged | grep ^cimageipi | awk -F"OS=" '{print $1}'| awk -F"\t" '{print $NF}' | cut -c1-50 | sed -e s/^\-/_/g > tmp.name_desc
		paste -d " " tmp.name_gene tmp.name_desc > tmp.name
	    else
		cat $p2.tagged | grep ^cimageipi | awk -F"[" '{print $1}'| awk -F"\t" '{print $NF}' | awk '{print $1}' > tmp.name_gene
		cat $p2.tagged | grep ^cimageipi | awk -F"[" '{print $1}'| awk -F"\t" '{print $NF}' | cut -d" " -f2- | cut -c1-50 | sed -e s/^\-/_/g > tmp.name_desc
		paste -d " " tmp.name_gene tmp.name_desc > tmp.name
	    fi
	fi
	paste tmp.ipi tmp.name >> tmp.ipi_name
    done
done

echo Creating input files for xcms
## create ipi number to protein name map
echo "name" > ipi_name.table
cat tmp.ipi_name | tr -d '\r' | sort | uniq | sed -e s/\'//g | sed -e s/\"//g | sed -e s/\;/\ /g >> ipi_name.table

## table with scans from different dataset
scanfiles=""
keys=$(cat tmp.key_scan | awk '{print $NF}' | sort | uniq )
for p11 in $(echo $mzxml)
do
    p1=$(echo $p11 | sed 's/_HL$//g')
    echo $p1 > $p1.tmp.scan
    for key in $keys
    do
	match=$(cat tmp.key_scan | grep -F $key | grep "^$p1 " | wc -l)
	#echo $key $match
	if [ "$match" != "0" ]; then
	    cat tmp.key_scan | grep -F $key | grep "^$p1 " | sort -k 4 -rn | head -1 | awk '{print $2}' >> $p1.tmp.scan
	else
	    echo "0" >> $p1.tmp.scan
	fi
    done

    scanfiles="$scanfiles $p1.tmp.scan"
done

echo "key mass" > tmp.seq_mass
for key in $keys
do
    cat tmp.key_scan | grep -F $key | head -1 | awk '{print $NF, $3}' >> tmp.seq_mass
done
paste -d " " tmp.seq_mass $scanfiles > cross_scan.table

## table with all ms2 scans
echo "key run scan HL" > all_scan.table
cat tmp.key_scan | awk '{print $NF, $1, $2, $(NF-1)}' >> all_scan.table

rm -rf tmp.ipi tmp.name* tmp.ipi_name tmp.key_scan tmp.seq_mass  *.tmp.scan

# filter by my_list.txt
mylist=$(\ls my_list.txt 2> /dev/null | wc -l)

if [ "$mylist" -eq 1 ];
then
    echo User provides a customized list in my_list.txt -- filtering...
    mv my_list.txt my_list.txt.original
    tr '\r' '\n' < my_list.txt.original > my_list.txt
    head -1 cross_scan.table > cross_scan.table.tmp
    for ml in $(cat my_list.txt);
    do
	cat cross_scan.table | grep ^$ml >> cross_scan.table.tmp
    done
    mv cross_scan.table cross_scan.table.all
    mv cross_scan.table.tmp cross_scan.table
fi

echo Running xcms to extract chromatographic peaks
echo "R --vanilla --args $param $mzxml < $CIMAGE_PATH/R/findMs1AcrossSetsFromDTASelect.R"
R --vanilla --args $param $mzxml < $CIMAGE_PATH/R/findMs1AcrossSetsFromDTASelect.R > findMs1AcrossSetsFromDTASelect.Rout

# echo Generating graphs
# cd output
# for p in $(\ls *.pdf | sed 's/\.pdf//g')
# do
#     ##ps2pdf $p.ps
#     if [ "$no_png" == "no_png" ]; then
# 	echo no png graphic conversion as requested
#     else
# 	##pdftops $p.pdf $p.ps
# 	mkdir -p PNG
# 	npages=$(cat $p.to_excel.txt | wc -l)
# 	nblock=$(($npages/500))
# 	echo $p.pdf has $npages pages
# 	##if [ $# -lt 3 ]; then
# 	##    rotate="-rotate 90"
# 	##else
# 	##    rotate=""
# 	##fi
# 	nb=0;
# 	while [ $nb -le $nblock ]; do
# 	    ns=$(($nb*500+1))
# 	    ne=$((($nb+1)*500))
# 	    if [ $ne -ge $npages ]; then ne=$(($npages-1)); fi
# 	    mkdir -p PNG/$nb
# 	    echo converting pages $ns-$ne
# 	   ## psselect -q -p$ns-$ne $p.ps ./PNG/$nb/$p.$nb.ps
# 	    pdftk A=$p.pdf cat A$ns-$ne output ./PNG/$nb/$p.$nb.pdf
# 	    cd PNG/$nb
# 	    ##convert  $p.$nb.ps $p.$nb.png
# 	    ##rm -rf $p.$nb.ps
# 	    pdftk $p.$nb.pdf burst
# 	    for pp in $(\ls pg_*.pdf | sed 's/\.pdf//g')
# 	    do
# 		pn=$(echo $pp | cut -f2 -d\_)
# 		pn2=$(echo $pn | awk '{print $1-1}')
# 		convert $pp.pdf $p.$nb\_$pn2.png 2> /dev/null
# 	    done
# 	    rm -rf *.pdf doc_data.txt
# 	    cd ../../
# 	    nb=$(($nb+1))
# 	done
# 	##rm -rf $p.ps
#     fi
#     echo done with $p.pdf
# done
# cd ..

mkdir -p output/TEXT
cp ipi_name.table output/TEXT
cp cross_scan.table output/TEXT
cp all_scan.table output/TEXT

echo Finished
exit 0
