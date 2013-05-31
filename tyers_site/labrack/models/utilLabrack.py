import os

from labrack.models.models import DnaComponent
  
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
            