from django.db import connections
from django.core.exceptions import ObjectDoesNotExist
from django.db.utils import ConnectionDoesNotExist
from models import PlasmidSample,Source, DnaComponent,DnaComponentType
from datetime import datetime
import MySQLdb
import sys
from django.db.models import Max

import re

def setup_cursor():
    try:
        cursor = connections['legacy'].cursor()
    except ConnectionDoesNotExist:
        print "Legacy database is not configured"
        return None
    
def decode_heuristically(string, enc = None, denc = sys.getdefaultencoding()):
    """
    Try to interpret 'string' using several possible encodings.
    @input : string, encode type.
    @output: a list [decoded_string, flag_decoded, encoding]
    """
    if isinstance(string, unicode): return string, 0, "utf-8"
    try:
        new_string = unicode(string, "ascii")
        return string, 0, "ascii"
    except UnicodeError:
        encodings = ["utf-8","iso-8859-1","cp1252","iso-8859-15"]

        if denc != "ascii": encodings.insert(0, denc)

        if enc: encodings.insert(0, enc)

        for enc in encodings:
            if (enc in ("iso-8859-15", "iso-8859-1") and
                re.search(r"[\x80-\x9f]", string) is not None):
                continue

            if (enc in ("iso-8859-1", "cp1252") and
                re.search(r"[\xa4\xa6\xa8\xb4\xb8\xbc-\xbe]", string)\
                is not None):
                continue

            try:
                new_string = unicode(string, enc)
            except UnicodeError:
                pass
            else:
                if new_string.encode(enc) == string:
                    return new_string, 0, enc

        # If unable to decode,doing force decoding i.e.neglecting those chars.
        output = [(unicode(string, enc, "ignore"), enc) for enc in encodings]
        output = [(len(new_string[0]), new_string) for new_string in output]
        output.sort()
        new_string, enc = output[-1][1]
        #return new_string, 1, enc
        return new_string
    
    
def to_unicode_or_bust(obj, encoding='utf-8'):
    if isinstance(obj, basestring):
        if not isinstance(obj, unicode):
            obj = unicode(obj, encoding)
    return obj

def whatisthis(s):
    if isinstance(s, str):
        print "ordinary string"
    elif isinstance(s, unicode):
        print "unicode string"
    else:
        print "not a string"
        
def import_plasmid():
    db = MySQLdb.connect(host="localhost", user="root", passwd="root",db="LIMS")
    cursor = db.cursor()
    
    
    if cursor is None:
        return
    ## it's important selecting the id field, so that we can keep the publisher - book relationship
    sql = """SELECT mtid, filing_date,source,plasmid_name,plasmid_marker, \
        plasmid_insert, plasmid_type, parental_vector, host_genotype, \
        attached_sheet, base_vector , synonyms, comments, disabled from plasmids"""
    cursor.execute(sql)
   
    DnaComponent.objects.all().delete()
    markerList = []
    for row in cursor.fetchall():
  
        if (row[0] <> '00000000'):
            try:
                try:                
                    mtid_value = decode_heuristically(row[0])
                    mtid_value = 'MTID'+mtid_value[0]
                    mtids = DnaComponent.objects.get(displayId=str(mtid_value))
                except:
                    sourcename = decode_heuristically(row[2])
                    source_old_list_val = getSource(sourcename[0]) 
                    
                    datevalue = decode_heuristically(row[1])
                    fillingDate_val = getDate(datevalue[0])
                    
                    markerValue = decode_heuristically(row[4])
                    markerValue = markerValue[0]
                    if ('cannot' not in markerValue):
                        markerValueList = []
                        if (markerValue.find(',')<>-1):
                            markerValueList = markerValue.split(",")
                        else:
                            if (markerValue.find(' ')<>-1):
                                markerValueList = markerValue.split(" ")
                            else:
                                if (markerValue.find('And')<>-1):
                                    markerValueList = markerValue.split("And")
                                else:
                                    if (markerValue.find('-')<>-1):
                                        markerValueList = markerValue.split("-")
                        
                            #if l not in markerList:
                                #markerList.append(l.replace(" ",""))
                    ##marker = getMarker(markerValue[0])
                    #plasmid = Plasmid(mtid=row[0], fillingDate=fillingDate_val, source_old_list=source_old_list_val, plasmidName=str(row[3]) \
                                      #, plasmidMarker=str(row[4]), plasmidInsert=str(row[5]), plasmidType=str(row[6]), parentalVector=str(row[7])\
                                      #, hostGenotype=str(row[8]), attachedSheet=str(row[9]), baseVector=str(row[10]), synonyms=str(row[11])\
                                      #, comments=str(row[12]))
                    testName = decode_heuristically(row[3])
                    testName = testName[0]
                    plasmid = DnaComponent(displayId=row[0], registrationDate=fillingDate_val, name=testName)
                    plasmid.save()
                    for l in markerValueList:
                        ##mark = 
                        mark = getMarker(l.upper())
                        plasmid.marker.add(mark)
                        
                    plasmid.save()
                    print str(row[0])
            except  Exception,e:
                #print "error :"+row[0]+" "+row[1]+" "+row[2]+" "+row[3]+" "+row[4]+" "+\
                      #row[5]+" "+row[6]+" "+row[7]+" "+row[8]+" "+row[9]+" "+\
                      #row[10]+" "+row[11]+" "+row[12]
                print "error :"+row[0]
                print str(e)
                pass
                 
    for m in markerList:
        print '---'
        print m
    cursor.close()
    db.close()
    
def getSource(value):
    if value<>'':
        try:
            source_old_list_val_value = decode_heuristically(value)
            source_old_list_val = Source.objects.get(name=source_old_list_val_value)
        except Source.DoesNotExist:
            source_old_list_val = Source(name = source_old_list_val_value)
            source_old_list_val.save()
    else:
        source_old_list_val = None    
    return source_old_list_val
def getDate(value):
    if (value.strip()!=''):
        if ("/" in value):
            d = datetime.strptime(value, '%d/%m/%Y')
            fillingDate_val = d.strftime('%Y-%m-%d')
        if ("." in value):
            d = datetime.strptime(value, '%d.%m.%Y')
            fillingDate_val = d.strftime('%Y-%m-%d')  
    else:
        d = datetime.strptime('01.01.1900', '%d.%m.%Y')
        fillingDate_val = d.strftime('%Y-%m-%d')      
    return fillingDate_val

def getMarkerType(valueType):
    if valueType<>'':
        try: 
            marker_value = valueType
            markerType = DnaComponentType.objects.get(name=marker_value)
        except DnaComponentType.DoesNotExist:                  
            markerTypeParent = DnaComponentType.objects.get(name = 'Marker')
            markerType = DnaComponentType(name=marker_value)
            markerType.save()
            markerType.subTypeOf=markerTypeParent
            markerType.save()
    else:
        markerType = None    
    return markerType

def getMarker(value):
    if value<>'':
        try:
            marker_value = decode_heuristically(value)
            marker_value = marker_value[0]
            marker_value_ws = marker_value
            marker_value_ws =marker_value_ws.replace(" ", "")                
            marker_value = DnaComponent.objects.get(displayId=marker_value_ws,componentSubType__subTypeOf__name='Marker')        
        except DnaComponent.DoesNotExist:                  
            dnaType = getMarkerType(value)
            marker_value = DnaComponent(displayId=marker_value_ws, name=marker_value,status='')
            marker_value.save()            
            marker_value.componentSubType.add(dnaType)
            marker_value.save()
    else:
        marker_value = None    
    return marker_value

def import_yeast():
    cursor = setup_cursor()
    if cursor is None:
        return
    sql = """SELECT title, company_id FROM books"""
    cursor.execute(sql)
    for row in cursor.fetchall():
        try:
            publisher = models.Publisher.objects.get(id=row[1])
        except ObjectDoesNotExist:
            print "Publisher not found with id %s" % row[1]
            continue
        else:
            book = models.Book(title=row[0], publisher=publisher)
            book.save()

def importMySQL():
    import_plasmid()
    print 'he'
