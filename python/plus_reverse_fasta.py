#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os, sys

fn = sys.argv[1]
#key_lst = []

handle = open( fn, "rU" )
full_handle = open( "full_"+fn, "w" )
for record in SeqIO.parse(handle, "fasta"):
    #if record.id in key_lst: continue
    #key_lst.append(record.id)
    SeqIO.write( record, full_handle, "fasta" )
    bw_rec = record
    bw_rec.id = "Reverse_" + record.id
    bw_rec.description = " ".join(record.description.split()[1:])
    bw_rec.seq = bw_rec.seq[::-1]
    SeqIO.write( bw_rec, full_handle, "fasta" )

handle.close()
