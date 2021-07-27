rm -rf *.pepindex
philosopher workspace --init
philosopher database --annotate target.fasta --prefix Reverse_
philosopher peptideprophet --nonparam --expectscore --decoyprobs --masswidth 1000.0 --clevel 2 --decoy Reverse_ --database target.fasta --combine $*
philosopher proteinprophet --maxppmdiff 2000000 --output combined interact.pep.xml
philosopher filter --sequential --razor --mapmods --prot 0.01 --tag Reverse_ --pepxml . --protxml combined.prot.xml
philosopher report --mzid
philosopher workspace --clean
