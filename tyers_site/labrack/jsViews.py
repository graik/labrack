from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from django.core import serializers

from django.utils import simplejson
from django.http import HttpResponse


import os
import labrack.models as M
import utilLabrack

from tyers_site import settings



def searchDnaParts(request, sequence_text,displayIdDnaComponent):
    coo = sequence_text
    rs = sequence_text.split('__')
    sequence_text = rs[0]
    sequence_vector = rs[1]
    f = len(sequence_vector)
    r = len(sequence_text)
    seq_exceptVector = sequence_text.replace(sequence_vector, '')
    s = len(seq_exceptVector)

    message = {"list_dnas": "", "extra_values": "","parttypes_values": "","optimizedfor_values": "", \
               "reverse_list_dnas": "","reverse_extra_values": ""}
    if request.is_ajax():
        # calculate the reverse complement sequence and duplicate sequence for better Vector matching
        sequence_textDuplicate = sequence_text + sequence_text
        my_seqDuplicate = Seq(sequence_textDuplicate, IUPAC.unambiguous_dna)
        revseqDuplicate = my_seqDuplicate.reverse_complement()                

        my_seq = Seq(sequence_text, IUPAC.unambiguous_dna)
        revseq = my_seq.reverse_complement()

        my_seq_exceptVect = Seq(seq_exceptVector, IUPAC.unambiguous_dna)
        revseq_exceptVect = my_seq_exceptVect.reverse_complement()

        message['list_dnas'] = getVectorBySequence(sequence_textDuplicate,'+',displayIdDnaComponent)
        message['reverse_list_dnas'] = getVectorBySequence(revseqDuplicate,'-',displayIdDnaComponent)

        #part for retriving potential Inserts
        message['extra_values'] = getInsertDBAnnotationBySequence(seq_exceptVector,'+',displayIdDnaComponent)
        message['reverse_extra_values'] = getInsertDBAnnotationBySequence(str(revseq_exceptVect),'-',displayIdDnaComponent)

        #paret retrieving all partTypes
        dnapartstypesAll = M.DnaComponentType.objects.all()
        json_dnaparttype = ''
        for dnaparttype in dnapartstypesAll :
            id = dnaparttype.id
            name = dnaparttype.name
            if json_dnaparttype == '':
                json_dnaparttype = '{ "id":"'+str(id)+'","name":"'+name+'"}'
            else:
                json_dnaparttype = json_dnaparttype+',{ "id":"'+str(id)+'","name":"'+name+'"}'
        json_dnaparttype = '['+json_dnaparttype+']'
        message['parttypes_values'] = json_dnaparttype

        #paret retrieving all optimized for
        chassisOptimizedAll = M.Chassis.objects.all()
        json_chassisOptimizedAll = ''
        json_chassisOptimizedAll = '{ "id":"","name":""}'
        for chas in chassisOptimizedAll :
            id = chas.id
            name = chas.name
            displayId = chas.displayId
            if json_chassisOptimizedAll == '':
                json_chassisOptimizedAll = '{ "id":"'+str(id)+'","name":"'+name+'","displayId":"'+displayId+'"}'
            else:
                json_chassisOptimizedAll = json_chassisOptimizedAll+',{ "id":"'+str(id)+'","name":"'+name+ \
                    '","displayId":"'+displayId+'"}'
        json_chassisOptimizedAll = '['+json_chassisOptimizedAll+']'
        message['optimizedfor_values'] = json_chassisOptimizedAll                
    else:
        message = "None"
    json = simplejson.dumps(message)
    return HttpResponse(json, mimetype='application/json')


def getDnaInfo(request, text_value):
    message = {"list_dnas": "", "extra_values": ""}
    try:
        if request.is_ajax():
            dnapart = M.DnaComponent.objects.get(displayId=text_value)
            json_object = ''
            id = dnapart.id
            name = dnapart.name
            displayid = dnapart.displayId
            sequence = dnapart.sequence
            description = dnapart.description
            json_object = '[{"id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+'","displayid":"'+displayid+ \
                '","description":"'+description+'"}]'
            message['list_dnas'] = json_object
            message['extra_values'] = request.method
        else:
            message = "You're the lying type, I can just tell."
        json = simplejson.dumps(message)
        return HttpResponse(json, mimetype='application/json')
    except:
        message = 'error'
        json = simplejson.dumps(message)
        return HttpResponse(json, mimetype='application/json')                




def saveUpload( uploaded, filename, raw_data ):
    '''
    raw_data: if True, uploaded is an HttpRequest object with the file being
              the raw post data
              if False, uploaded has been submitted via the basic form
              submission and is a regular Django UploadedFile in request.FILES
    '''
    #filename = settings.MEDIA_ROOT + '/documents/GenBank/'+ filename
    filename = settings.MEDIA_ROOT + '/'+ filename
    try:
        from io import FileIO, BufferedWriter
        with BufferedWriter( FileIO( filename, "wb" ) ) as dest:
            # if the "advanced" upload, read directly from the HTTP request
            # with the Django 1.3 functionality
            if raw_data:
                foo = uploaded.read( 1024 )
                while foo:
                    dest.write( foo )
                    foo = uploaded.read( 1024 )
            # if not raw, it was a form upload so read in the normal Django chunks fashion
            else:
                for c in uploaded.chunks( ):
                    dest.write( c )
            # got through saving the upload, report success
            return True
    except IOError:
        # could not open the file most likely
        return HttpResponseBadRequest( "Could not open the file" )
    return False

def fileUpload( request ):
    if request.method == "POST":   
        if request.is_ajax( ):
            # the file is stored raw in the request
            upload = request
            is_raw = True
            # AJAX Upload will pass the filename in the querystring if it is the "advanced" ajax upload
            try:
                filename = request.GET[ 'qqfile' ]
            except KeyError:
                return HttpResponseBadRequest( "AJAX request not valid" )
        # not an ajax upload, so it was the "basic" iframe version with submission via form
        else:
            is_raw = False
            if len( request.FILES ) == 1:
                # FILES is a dictionary in Django but Ajax Upload gives the uploaded file an
                # ID based on a random number, so it cannot be guessed here in the code.
                # Rather than editing Ajax Upload to pass the ID in the querystring,
                # observer that each upload is a separate request,
                # so FILES should only have one entry.
                # Thus, we can just grab the first (and only) value in the dict.
                upload = request.FILES.values( )[ 0 ]
            else:
                raise Http404( "Bad Upload" )
            filename = upload.name

        # save the file
        success = saveUpload( upload, filename, is_raw )
        genbankjson = retrieveGenBankInfo(filename)
        # let Ajax Upload know whether we saved it or not
        import json
        ret_json = { 'success': success, 'test':genbankjson,}
        return HttpResponse(json.dumps(ret_json))




def getGBDataFromFile(request, filePath):
    success = ''
    genbankjson = retrieveGenBankInfo(filePath)
        # let Ajax Upload know whether we saved it or not
    import json
    ret_json = { 'success': success, 'test':genbankjson,}
    return HttpResponse(json.dumps(ret_json))       




def getVectorBySequence(sequence_text,strand,displayIdDnaComponent):#local function

    dnapartsVectorAll = M.DnaComponent.objects.filter(componentType__name='Vector Backbone')
    dnVectorId = -999 
    try:
        dnaComp = M.DnaComponent.objects.get(displayId =displayIdDnaComponent)
        qSet = dnaComp.related_annotations()
        for dn in qSet:
            if dn.componentAnnotated in dnapartsVectorAll:
                dnVectorId = dn.componentAnnotated.id                
    except M.DnaComponent.DoesNotExist:
        dnVectorId = -999       


    dnaparts = [dnapart for dnapart in M.DnaComponent.objects.all() \
                if (dnapart in dnapartsVectorAll and dnapart.sequence is not None \
                    and dnapart.sequence != "" and dnapart.sequence in sequence_text)]

    name = ''
    data = serializers.serialize('json', dnaparts)
    json_object = ''
    isRelated = 'False'
    for dnapart in dnaparts :
        id = dnapart.id
        if (dnapart.id==dnVectorId):
            isRelated = 'True'
        else:
            isRelated = 'False'
        name = dnapart.name
        displayid = dnapart.displayId
        sequence = dnapart.sequence
        description = dnapart.description
        if json_object == '':
            json_object = '{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+ \
                '","displayid":"'+displayid+'","description":"'+description+'","strand":"'+strand+ \
                '","isRelated":"'+isRelated+'"}'
        else:
            json_object = json_object+',{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+ \
                '","displayid":"'+displayid+'","description":"'+description+'","strand":"'+strand+ \
                '","isRelated":"'+isRelated+'"}'
    json_object = '['+json_object+']'
    return json_object


def getAnnotToBeDeleted(request, jsonAmmpt,displayIdDnaComponent):
    import django.utils.simplejson as json
    try:
        dnacomp = M.DnaComponent.objects.get(displayId = displayIdDnaComponent)
    except M.DnaComponent.DoesNotExist:
        dnacomp = None
        message = "None"
        json = simplejson.dumps(message)
        return HttpResponse(json, mimetype='application/json')                

    try:
        data2 = json.loads('['+jsonAmmpt+']')
        sumDna = ""
        missingAnnot = ''
        selectedAnnotFromDB = json.loads(data2[0]["selected_annot"][0]["db_checked"])
        if (selectedAnnotFromDB):                     
            for id in selectedAnnotFromDB:
                sumDna = sumDna+','+id["id"]       
        try:
            dnacomp = M.DnaComponent.objects.get(displayId = displayIdDnaComponent)
            for qannot in M.DnaSequenceAnnotation.objects.filter( subComponent=dnacomp.id):
                test = qannot.componentAnnotated.displayId
                if( not utilLabrack.isDnaTypeVector(qannot.componentAnnotated)):
                    if (sumDna.find(test)==-1):
                        if (missingAnnot ==''):
                            missingAnnot = test
                        else:
                            missingAnnot = missingAnnot + "," + test
        except M.DnaComponent.DoesNotExist:
            pass
        json = simplejson.dumps(missingAnnot)
        return HttpResponse(json, mimetype='application/json') 
    except ValueError: ## includes simplejson.decoder.JSONDecodeError
        dnacomp = None
        message = "None"
        json = simplejson.dumps(message)
        return HttpResponse(json, mimetype='application/json')                 



def getInsertDBAnnotationBySequence(sequence_text,strand,displayIdDnaComponent):#local function
    lenSimple = len(sequence_text)        
    sequence_text = sequence_text + sequence_text

    dnapartsVectorAll  = M.DnaComponent.objects.filter(componentType__name='Vector Backbone')
    dnapartsInsertsAll = M.DnaComponent.objects.filter(componentType__name='Vector Backbone')
    # removed because Vector selection was removed for now, therefore add Vector selection within the annotations
    #dnaparts = [dnapart for dnapart in DnaComponent.objects.all() if (dnapart not in dnapartsVectorAll and dnapart.sequence is not None and dnapart.sequence != "" and dnapart.sequence in sequence_text)]
    dnaparts = [dnapart for dnapart in M.DnaComponent.objects.all()  \
                if (dnapart.sequence is not None and dnapart.sequence != "" and \
                    dnapart.sequence in sequence_text and dnapart.displayId!=displayIdDnaComponent)]
    name = ''
    data = serializers.serialize('json', dnaparts)
    json_Insertobject = ''
    for dnapart in dnaparts :
        id = dnapart.id

        name = dnapart.name
        displayid = dnapart.displayId
        sequence = dnapart.sequence
        fisrtPeace = len(sequence)
        secondpeace = len(sequence_text)
        rslt = fisrtPeace/float(secondpeace)
        firstposition = sequence_text.find(sequence)
        # check if the sequence was found over circular dna sequence by cmopare the lentgh of simpleSequence and not duplicate
        lastposition = firstposition+len(sequence)
        if (firstposition+len(sequence)>lenSimple):
            lastposition = fisrtPeace - (lenSimple - firstposition)


        # check if the annotation are already related
        isRelated = False
        isRelated = M.DnaSequenceAnnotation.isRelated(displayid,displayIdDnaComponent)
        
        if (strand=='-'):
            firstPos = firstposition
            firstposition = secondpeace - lastposition
            lastposition = secondpeace - firstPos
        # remove case where sequence from DB is the same as the requested annotation to avoid annotating fragment by it self
        if (sequence!=sequence_text):
            coverage = str(firstposition)+'-'+str(lastposition)
            description = dnapart.description
            if (dnapart.optimizedFor!=None):
                optimizedFor_name = dnapart.optimizedFor.name
            else:
                optimizedFor_name = ''
            if (dnapart.componentType!=None):
                componentType_name  = ''
                for cpType in dnapart.componentType.all():
                    if (componentType_name==''):
                        componentType_name = componentType_name+cpType.name
                    else:
                        componentType_name = componentType_name+','+cpType.name
            else:
                componentType_name = ''                        
            if json_Insertobject == '':
                json_Insertobject = '{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+ \
                    '","displayid":"'+displayid+'","description":"'+description+'","coverage":"'+coverage+ \
                    '","optimizedFor_name":"'+optimizedFor_name+'","componentType_name":"'+componentType_name+ \
                    '","strand":"'+strand+'","isRelated":"'+str(isRelated)+'"}' 
            else:
                json_Insertobject = json_Insertobject+',{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+\
                    '","displayid":"'+displayid+'","description":"'+description+'","coverage":"'+coverage+\
                    '","optimizedFor_name":"'+optimizedFor_name+'","componentType_name":"'+componentType_name+ \
                    '","strand":"'+strand+'","isRelated":"'+str(isRelated)+'"}' 
    json_Insertobject = '['+json_Insertobject+']'
    return json_Insertobject



def retrieveGenBankInfo(filename):
    #gb_file = settings.MEDIA_ROOT+"documents/GenBank/"+os.path.normpath(filename)
    gb_file = settings.MEDIA_ROOT+"/"+os.path.normpath(filename)
    gb_features = ""
    dispId = 1
    isParsingDone = False
    json_Insertobject = '';
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        if (not isParsingDone):
            for ind in xrange(len(gb_record.features)) :
                nameType = repr(gb_record.features[ind].type).replace("'", "")
                isParsingDone = True
                #gb_features += '\n'+ repr(gb_record.features[ind].type) + " Location start : "+ repr(gb_record.features[ind].location._start.position) + " Location end : "+ repr(gb_record.features[ind].location._end.position)
                strandValue = repr(gb_record.features[ind].strand)
                startPos = repr(gb_record.features[ind].location._start.position)
                endPos = repr(gb_record.features[ind].location._end.position)
                pos = str(startPos)+'-'+str(endPos)
                strandValue = repr(gb_record.features[ind].strand)
                origFullSequence = gb_record.seq.tostring()
                featureSeq = gb_record.features[ind].extract(gb_record).seq.tostring()
                coverage = pos
                description = ''
                label = repr(gb_record.features[ind].qualifiers.get('label')).replace("['","").replace("']","").replace("\\"," ")
                if (label == 'None'):
                    label = repr(gb_record.features[ind].qualifiers.get('gene')).replace("['","").replace("']","").replace("\\"," ")

                if json_Insertobject == '':
                    json_Insertobject = '{ "id":"'+str(label)+'","name":"'+nameType+'","sequence":"'+origFullSequence+ \
                        '","displayid":"'+label+'","description":"'+description+'","coverage":"'+coverage+ \
                        '","startPos":"'+startPos+'","endPos":"'+endPos+'","strandValue":"'+strandValue+'","featureSeq":"'+ \
                        featureSeq+'"}' 
                else:
                    json_Insertobject = json_Insertobject+',{ "id":"'+str(label)+'","name":"'+nameType+'","sequence":"'+ \
                        origFullSequence+'","displayid":"'+label+'","description":"'+description+'","coverage":"'+coverage+ \
                        '","startPos":"'+startPos+'","endPos":"'+endPos+'","strandValue":"'+strandValue+ \
                        '","featureSeq":"'+featureSeq+'"}' 
    json_Insertobject = '['+json_Insertobject+']' 
    return json_Insertobject