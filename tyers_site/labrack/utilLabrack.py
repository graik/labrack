import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
            
def getStrand(fullSeq,seq):

    first = fullSeq.find(seq)
    if first >= 0:
        return '1'

    my_seq = Seq(seq, IUPAC.unambiguous_dna)
    revseq = my_seq.reverse_complement()     ## import this module
    first = fullSeq.find(str(revseq))
    if first >= 0:
        return '-1'
    else:
        return '?'
            