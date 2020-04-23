#!/usr/bin/perl -w

die "Usage: $0 text_table [ratio.png] [run_dirs] \n" if (@ARGV <1 );

$intable="$ARGV[0].txt";
$cwd=$ARGV[1];
$dset="$cwd/$intable";
open(INFILE,$intable) || die "cannot open $intable: $!\n";
@txt = <INFILE>;
close(INFILE);
$nrow = @txt;
@header = split(/\t/,$txt[0]);
$ncol = @header;
@bgcolormap=@header;
$nset=0;
for ($i=0; $i<$ncol; ++$i) {
    if ($header[$i] =~ /^mr\./) {
	$bgcolormap[$i]="bgcolor=\"#DCDCDC\"";
	$nset++;
    } else {
	$bgcolormap[$i]=""
    }
}

$hostname = "162.105.22.250";

$outhtml="$ARGV[0].html";
##$outhtml =~ s/\.txt$/\.html/g;

@cwdstr = split("\/",$cwd);
$str1 = pop(@cwdstr);
$str2 = pop(@cwdstr);
$htmltitle = join(":",$str2, $str1);
open(OUTFILE, ">$outhtml") || die "cannot write to $outhtml: $!\n";
##print OUTFILE "Content-type: text/html\n\n";
print OUTFILE <<ENDOFHEADER;
<HTML><HEAD>
<STYLE TYPE="text/css">
*{
    font-family: arial, helvetica, myriad;
    font-size: 13px;
}
table#sample TD {
    text-align: center;
}
</STYLE>
<TITLE>
$htmltitle
</TITLE>
</HEAD>
<BODY>
ENDOFHEADER

print OUTFILE "Data Location: $cwd (<A HREF=\"./\">Go To</A>)<br><br>\n";
print OUTFILE "<A HREF=\"$ARGV[0]\.txt\">Tab-delimited Text Format</A><br><br>\n";
print OUTFILE "<A HREF=\"$ARGV[0]\.to_excel.txt\">Ratio Comparsion Table(to Excel)</A><br><br>\n";

print OUTFILE "<A HREF=\"http://$hostname/cgi-bin/chuquest/batch_annotate.pl?dset=$cwd/$intable\">Gabe's Uniprot Batch Annotation</A><BR><BR>\n";
print OUTFILE "<A HREF=\"$ARGV[0]\.png\">Ratio Plot</A><BR><BR>\n";
if( -e "$ARGV[0].vennDiagram.png" ) {
    print OUTFILE "<A HREF=\"$ARGV[0]\.vennDiagram\.png\">Venn Diagram</A><BR><BR>\n";
}
if(@ARGV>2) {
    for ( $i=2; $i<@ARGV; $i++ ) {
	$j = $i-1;
	print OUTFILE "run $j: <A HREF=\"$ARGV[$i]\">$ARGV[$i]</A><BR><BR>";
  }
}
print OUTFILE "<b>mr.set columns:</b><BR>";
print OUTFILE "[bold] -- median value of measured ratios in this group<BR>";
print OUTFILE "<b>sd.set columns:</b><BR>";
print OUTFILE "[bold] -- standard deviation of measured ratios in this group<BR>";
print OUTFILE "[plain] -- for each labeled peptide, number of peaks with triggered ms2 / number of candidate peaks for selection<BR>";
print OUTFILE "<b>ratio calculation:</b><BR>";
print OUTFILE "peak pairs with R2<0.8 are excluded<BR>";
print OUTFILE "labeled peptides with invalid 1:10 ratios are excluded<BR>";
print OUTFILE "<b>peptide grouping:</b><BR>";
print OUTFILE "multiply labeled peptides are separated from singly labeled ones<BR>";
print OUTFILE "singly labeled peptides with multiple alternative sites are grouped together<BR><BR>";

print OUTFILE "<TABLE id=\"sample\" border=2 frame=\"border\" rules=\"groups\" summary=$intable>\n";
##print OUTFILE "<CAPTION> <b>$title</b> </CAPTION>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=4>\n";
#print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=$nset>\n";
print OUTFILE "<COLGROUP span=$nset>\n";
#print OUTFILE "<COLGROUP span=3>\n";
print OUTFILE "<COLGROUP span=1>\n";
#print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<TBODY>\n";
print OUTFILE "<TR>\n";
for ($i=0; $i<$ncol; ++$i) {
    #print OUTFILE "<TH align=\"center\" bgcolor=\"HoneyDew\">$header[$i]\n";
    print OUTFILE "<TH bgcolor=\"HoneyDew\">$header[$i]\n";
}
print OUTFILE "<TR>\n";
@line1 = split('\t', $txt[1]);
for ($i=0; $i<$ncol; ++$i) {
    if ( $line1[$i] eq ' ' ) {
	print OUTFILE "<TH bgcolor=\"HoneyDew\"> \n";
} else {
    print OUTFILE "<TH bgcolor=\"HoneyDew\"><A HREF=\"http://$hostname/cgi-bin/chuquest/cimage_sort_compare_table.pl?dset=$dset&colname=$header[$i]&ascending=True\" style=\"text-decoration:none\">^</A> <A HREF=\"http://$hostname/cgi-bin/chuquest/cimage_sort_compare_table.pl?dset=$dset&colname=$header[$i]&ascending=False\" style=\"text-decoration:none\">v</A>\n";
}
}
print OUTFILE "<TBODY>\n";
for ($i=1; $i<$nrow; ++$i) {
    print OUTFILE "<TR>";
    @line = split('\t',$txt[$i]);
    if ($line[0] =~ /\d+/ ) {
	$bold1="<b>";
	$bold2="</b>";
	$anchor1 = 1;
    } else {
	$bold1="";
	$bold2="";
	$anchor1 = 0;
    }
    for ($j=0; $j<$ncol; ++$j) {
	$_ = $line[$j];
	if ( /^=HYPERLINK/ ) {
	    /^=HYPERLINK\("(.*)","(\d+.\d+.\d+)"\)$/;
	    print OUTFILE "<TD $bgcolormap[$j]><A HREF=\"$1\">$2</A>";
	} elsif ( /^([a-zA-Z]+)$/ && $anchor1 ) {
	    print OUTFILE "<TD $bgcolormap[$j]> <A NAME=\"$1\"></A> $bold1 $1 $bold2";
	} elsif ( /(\w+)$/ && $j==1 ) {
	    if (/^(IPI\w+)/) {
		print OUTFILE "<TD $bgcolormap[$j]><A HREF=\"http://www.ebi.ac.uk/cgi-bin/dbfetch?db=IPI&id=$1&format=default\">$1</A>";
	    } elsif (/^(\w{6})$/) {
		print OUTFILE "<TD $bgcolormap[$j]><A HREF=\"http://www.uniprot.org/uniprot/$1\">$1</A>";
	    } else {
		print OUTFILE "<TD $bgcolormap[$j]> $bold1 $line[$j] $bold2";
	    }
	} elsif (/^(\d+\/\d+)/) {
	    print OUTFILE "<TD $bgcolormap[$j]> $bold1 $1 $bold2";
	} else {
	    print OUTFILE "<TD $bgcolormap[$j]> $bold1 $line[$j] $bold2";
	}
    }
    print OUTFILE "\n";
}
print OUTFILE '</TABLE>'."\n";
print OUTFILE '</BODY>'."\n";
print OUTFILE '</HTML>'."\n";
