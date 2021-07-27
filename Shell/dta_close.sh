rm *.pepindex
philosopher workspace --init
philosopher database --annotate target.fasta --prefix Reverse_
philosopher peptideprophet --nonparam --expectscore --decoyprobs --ppm --accmass --decoy Reverse_ --database target.fasta *.pepXML 
philosopher proteinprophet --maxppmdiff 2000000 --output combined *.pep.xml
philosopher filter --sequential --razor --mapmods --prot 0.01 --tag Reverse_ --pepxml . --protxml combined.prot.xml
philosopher report --mzid
philosopher workspace --clean
