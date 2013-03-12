import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from labrack.models.component import DnaComponent

            
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
            
            
def getNextAvailableDNAName(userInitial):
    displayIdDefault = DnaComponent.objects.filter(displayId__startswith=userInitial).order_by('-displayId')[:1]
    if not displayIdDefault:
        displayIdDefaultName = ''
    else:
        dis = displayIdDefault[0].displayId
        initial = dis[:2]
        numberSeq = int(dis[2:6])+1
        numSeqChar = "0000"+str(numberSeq)
        numSeqChar = numSeqChar[1:6]
        displayIdDefaultName = initial+str(numSeqChar)+'a'
    
        
    return displayIdDefaultName
            