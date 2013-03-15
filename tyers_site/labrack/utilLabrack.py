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
        fullSeq = fullSeq + fullSeq
        first = fullSeq.find(seq)
        if first >= 0:
            return '1'
    
        my_seq = Seq(seq, IUPAC.unambiguous_dna)
        revseq = my_seq.reverse_complement()     ## import this module
        first = fullSeq.find(str(revseq))
        if first >= 0:
            return '-1'
        else:
            return '+1'
            
def isDnaTypeVector(dnaComponent):
    dnapartsVectorAll = DnaComponent.objects.filter(componentType__name='Vector Backbone')
    if dnaComponent in dnapartsVectorAll:
        return True
    
    return False
    
            
def getNextAvailableDNAName(currUser):
    lName  = currUser.last_name
    fName = currUser.first_name
    initialUser = lName[:1]+fName[:1]
    if initialUser=='':
        initialUser = str(currUser)[:2]
        
    displayIdDefault = DnaComponent.objects.filter(displayId__startswith=initialUser).order_by('-displayId')[:1]
    
    if not displayIdDefault:
        displayIdDefaultName = initialUser+'0001a'
    else:
        dis = displayIdDefault[0].displayId
        initial = dis[:2]
        numberSeq = int(dis[2:6])+1
        numSeqChar = "0000"+str(numberSeq)+"0"
        numSeqChar = str(numSeqChar)[-5:-1]
        displayIdDefaultName = initial+str(numSeqChar)+'a'
    
        
    return displayIdDefaultName
            