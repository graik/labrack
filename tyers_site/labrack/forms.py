from django import forms
from django.core.exceptions import ValidationError
from django.contrib import messages
from labrack.models.sample import DnaSample
from labrack.models.sample import ChassisSample
from labrack.models.generalmodels import Document
from labrack.models.component import DnaComponent
from labrack.models.component import Component
from labrack.models.component import DnaComponentType
from labrack.models.component import Chassis
from labrack.models.generalmodels import DnaSequenceAnnotation
from labrack.models.generalmodels import ProteinSequenceAnnotation
from django.core.files import File
from django.core.exceptions import ValidationError
from django.forms.fields import DateField, ChoiceField, MultipleChoiceField
#from django.forms.widgets import RadioSelect, CheckboxSelectMultiple
from django.forms.extras.widgets import SelectDateWidget
from django.forms.widgets import *
import django.utils.simplejson as json
from tyers_site import settings
from django.contrib.messages.storage.fallback import FallbackStorage
import os
import re

import utilLabrack



from django.forms.forms import NON_FIELD_ERRORS



class DnaSampleForm(forms.ModelForm):



    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)

        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(DnaSampleForm, self).__init__(*args, **kwargs)

    def clean(self):
        cleanedData = super(DnaSampleForm, self).clean()
        dataDnaConstruct = cleanedData.get("dnaConstruct")
        dataDnaDisplayID = cleanedData.get("DNA_Display_ID")
        dataDnaName = cleanedData.get("DNA_Name")
        dataDnaDescription = cleanedData.get("DNA_Description")
        che = cleanedData.get("inChassis")
        if dataDnaDisplayID=='':

            dataDnaDisplayID=None
        if (dataDnaConstruct==None and dataDnaDisplayID==None):
            raise forms.ValidationError("Either select an existing DNA Construct or fill the dna description to create a new one!")
        return cleanedData        


    def clean_DNA_Display_ID(self):
        data = self.cleaned_data['DNA_Display_ID']
        if (data!= '' and DnaComponent.objects.filter(displayId=data)):
            raise forms.ValidationError("This DNA Display ID:"+data+" Already exists!")
        return data

    def save(self, commit=True):
        instance = super(DnaSampleForm, self).save(commit=False)
        #instance.historyDescription = self.cleaned_data['genebank'] # etc
        
        if self.cleaned_data['dnaConstruct']==None:
            dna2db = DnaComponent(displayId=self.cleaned_data['DNA_Display_ID'],description=self.cleaned_data['DNA_Description'], name =self.cleaned_data['DNA_Name'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['DNA_Status'], genBankfile = self.cleaned_data['genebank'])
            dna2db.saveSequenceWithoutAnnotations()
            #dna2db.save()
            instance.dnaConstruct = dna2db
        else:
            instance.dnaConstruct = self.cleaned_data['dnaConstruct']
        if commit:
            instance.save()
        return instance 

    class Meta:
        model = DnaSample


class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a CSV File'
    )

class ChassisSampleForm(forms.ModelForm):

    ChassisDisplayID = forms.CharField(required=False,label="Display ID")
    ChassisName = forms.CharField(required=False,label="Name")  
    ChassisDescription = forms.CharField(widget=forms.TextInput(attrs={'size':'40'}),required=False,label="Description")

    historyDescription = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False,label="Description")
    class Meta:
        model = ChassisSample    
        
        
        




class DnaComponentForm(forms.ModelForm):


    htmlAttribute1 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800'}),required=False,label="")
    htmlAttribute2 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800','hidden':'true'}),required=False,label="")
    htmlAttribute3 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800','hidden':'true'}),required=False,label="")
    
    errorMessageForm = ''
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(DnaComponentForm, self).__init__(*args, **kwargs)
        
    
        
    def clean(self):
        errorMessage = ""
        errors = {}        
        
        jsonCancel = self.cleaned_data['htmlAttribute2']
        try:
            dataCancel = json.loads(jsonCancel)
           
            gb_checked_cancel = json.loads(dataCancel["selected_annot"][2]["gb_checked_cancel"])        
            if (gb_checked_cancel=='True'):
                errorMessage = "The user has cancel the saving procedure due to annotation deletions risks!"
                errors.setdefault('',[]).append(errorMessage)                
        except ValueError:
            ## if there is no valid json data yet, just pass
            pass
        
        # check if disp name is clean with pattern XXX1234X
        try:
            dispId = self.cleaned_data['displayId']
        except Exception, err:
            errorMessage = "ID field is mandatory!"
            errors.setdefault('',[]).append(errorMessage)                      
            raise forms.ValidationError("ID field is mandatory!")      
            commit=False
       
        dispId = self.cleaned_data['displayId']
        compType = self.cleaned_data['componentType']
        compTypeName = 'Others DNA Type'        
        for s in compType:
            if s.name =='Vector Backbone':
                compTypeName = 'Vector Backbone'
            if s.name =='Primer':
                compTypeName = 'Primer'                
        
            
        regExp = r'^[a-zA-Z][a-zA-Z]\d\d\d\d[a-zA-Z]?$'
        regExp2 = r'^[a-zA-Z][a-zA-Z]\d\d\d\d?$'
        if (compTypeName=='Vector Backbone'):
            regExp = r'^[vV][a-zA-Z]\d\d\d?$'
            regExp2 = r'^[vV][a-zA-Z]\d\d\d?$'
        if (compTypeName=='Primer'):
            regExp = r'^O[a-zA-Z][a-zA-Z]\d\d\d\d[a-zA-Z]?$'
            regExp2 = r'^O[a-zA-Z][a-zA-Z]\d\d\d\d?$'            
            
        if not(re.match(regExp, dispId)):
            if (compTypeName=='Vector Backbone'):
                errorMessage = "ID :"+dispId+" doesn't respect the naming convention for "+compTypeName
                errors.setdefault('',[]).append(errorMessage)
            else:
                if (compTypeName=='Primer'):
                    errorMessage = "ID :"+dispId+" doesn't respect the naming convention for "+compTypeName
                    errors.setdefault('',[]).append(errorMessage)
                else:
                    ## raise forms.ValidationError("ID :"+dispId+" doesn't respect the pattern XX1234X, 
                    ## example : sb0001A where sb for initial and A for version")
                    errorMessage = "ID :"+dispId+" doesn't respect the naming convention for "+compTypeName
                    errors.setdefault('',[]).append(errorMessage)                      
        
        
        cleanedData = super(DnaComponentForm, self).clean()
        
        ## retrieve vector GB information
        jsonFromWeb = self.cleaned_data['htmlAttribute1']
        fullSequence = self.cleaned_data['sequence']
     
        try:            
            ## retrieve annot GB information
            jsonFromWeb2 = self.cleaned_data['htmlAttribute2']
            data2 = json.loads(jsonFromWeb2)
            selectedAnnotFromGB = json.loads(data2["selected_annot"][1]["gb_checked"])
            
            if (self.instance.id!=None):        
                r = self.instance.number_related_ParentChildAnnotations()
                if (fullSequence<>self.instance.sequence and self.instance.sequence<>"" and r>0) :
                    
                    print ''
        except ValueError:
            selectedAnnotFromGB=[]
            self.errorMessageForm = 'no json data was submited and no annotation data was saved due to errors, contact administrator!' 
            print 'no GB json data was submited!'          
               
           
        # check if the named annotation exist already in the DB
        if (selectedAnnotFromGB):                
                
                for obj in selectedAnnotFromGB:
                    
                    idAnnot=obj["text2"]
                    annotType=obj["text5"]
                    regExp = r'^[a-zA-Z][a-zA-Z]\d\d\d\d[a-zA-Z]?$'
                    regExp2 = r'^[a-zA-Z][a-zA-Z]\d\d\d\d?$'
                    
                    if (annotType=='Primer'):
                        regExp = r'^O[a-zA-Z][a-zA-Z]\d\d\d\d[a-zA-Z]?$'
                        regExp2 = r'^O[a-zA-Z][a-zA-Z]\d\d\d\d?$'                            
                    
                    if (idAnnot.upper()<>'(AUTO)'):
                        if (not(re.match(regExp, idAnnot)) and not(re.match(regExp2, idAnnot))):
                            errorMessage = "ID Annotations Error: "+idAnnot+" doesn't respect the naming convention for "+annotType
                            errors.setdefault('',[]).append(errorMessage)                        
                        
                    if (DnaComponent.objects.filter(displayId=idAnnot)):
                        errorMessage = "ID Annotations Error: Please reload the GB file and choose another name for the Annotation "+idAnnot+" since it already exist in the DB!"
                        errors.setdefault('',[]).append(errorMessage)
                        
                    if (dispId==idAnnot):
                        errorMessage = "ID Annotations Error: the name for the Annotation "+idAnnot+" should be different from the created DNA Part!"
                        errors.setdefault('',[]).append(errorMessage) 
        if  (len(errors)>0):
            raise forms.ValidationError(errors)
        print '+ methods clean is done'
        return cleanedData 

    def get_success_message(self, cleaned_data):
        return self.success_message % dict(cleaned_data,
                                           calculated_field=self.object.calculated_field)
    
    def save(self, commit=True):
        print '+ methods done for dnaComponent'
        instance = super(DnaComponentForm, self).save(commit=False)
        fileName = self.cleaned_data['htmlAttribute3']
        
        if (fileName!=""):
            #file = open(settings.MEDIA_ROOT+"/documents/GenBank/"+fileName)
            file = open(settings.MEDIA_ROOT+"/"+fileName)
            instance.genBankfile = File(file)

        jsonFromWeb = self.cleaned_data['htmlAttribute1']
        
        fullSequence = self.cleaned_data['sequence']
        try:
            created_by = self.cleaned_data['created_by']
        except KeyError,NameError:
            created_by = instance.created_by
       
        jsonFromWeb2 = self.cleaned_data['htmlAttribute2']
        try:
            data2 = json.loads(jsonFromWeb2)
            selectedAnnotFromDB = json.loads(data2["selected_annot"][0]["db_checked"])
            selectedAnnotFromGB = json.loads(data2["selected_annot"][1]["gb_checked"])
        except ValueError:
            selectedAnnotFromDB=[]
            selectedAnnotFromGB=[]
            self.errorMessageForm = 'no json data was submited and no annotation data was saved due to errors, contact administrator!' 
            print 'no  json data was submited!'  

        isVectorBackbone = True
        if selectedAnnotFromDB and selectedAnnotFromGB:
            isVectorBackbone = False

        #get the parttypeVector
        if (not DnaComponentType.objects.filter(name='Vector Backbone')):
            subCtType = DnaComponentType(name = 'Vector Backbone')
            subCtType.save()
        subCtVectorType = DnaComponentType.objects.filter(name='Vector Backbone')

        
        instance.save()
        #delete all the entry for annotation for this component ID and recreate if needed
        DnaSequenceAnnotation.deleteRelated(instance.displayId)        
        id = DnaComponent.objects.get(displayId=instance.displayId)                
        coun = DnaSequenceAnnotation.objects.filter( subComponent=id)

        print '+ start selectedAnnotFromDB'
        if (selectedAnnotFromDB):
            try:                       
                for id in selectedAnnotFromDB:
                    dnaAnnot = DnaComponent.objects.get(displayId=id["id"])
                    coverage=id["text6"].replace('(-)','').replace('(+)','')
                    coverage = coverage.split('-')
                    first = int(coverage[0])
                    secon = int(coverage[1])                            
                    stra = utilLabrack.getStrand(fullSequence,dnaAnnot.sequence)
                    an2db = DnaSequenceAnnotation(uri ='', bioStart = first, bioEnd = secon, strand = stra, subComponent = instance, componentAnnotated = dnaAnnot)
                    an2db.save()

            except M.DnaComponent.DoesNotExist:
                commit=False
        print '+ end selectedAnnotFromDB'
        print '+ start selectedAnnotFromGB'
        if (selectedAnnotFromGB):                
            for obj in selectedAnnotFromGB:
                coverage=obj["text6"]
                id=obj["text2"]
                coverage = coverage.split('-')
                first = int(coverage[0])
                secon = int(coverage[1])
                #seq = fullSequence[first:secon]
                if (first>secon):
                    fir = fullSequence[first:]
                    sec = fullSequence[:secon]
                    seq = fir+sec
                else:
                    seq = fullSequence[first:secon]
                name=obj["text3"]
                descrp=obj["text4"]
                dnatype=obj["text5"]
                optimizedfor=obj["text7"].split(' - ')
                try:
                    chassisOptimizedFor = Chassis.objects.get(displayId=optimizedfor[0])
                except Chassis.DoesNotExist:
                    chassisOptimizedFor = None
                subCtVectorType = DnaComponentType.objects.filter(name__in=[dnatype,'Annotation'])	

                dnaAnnot = DnaComponent(displayId=id,description=descrp, name =name, sequence = seq,optimizedFor = chassisOptimizedFor, created_by = created_by)
                dnaAnnot.save()
                if (id.upper()=='(AUTO)' or id.upper()==''):
                    
                    if (dnatype=='Vector Backbone' or dnatype=='Vector'):
                        lastComponent = DnaComponent.objects.filter(displayId__regex='^Ve.*').order_by('-displayId')[:1]
                        leng = len(lastComponent)
                        if leng>0:
                            lastComponentDisplayID = int(lastComponent[0].displayId[2:])
                        else:
                            lastComponentDisplayID = 0                                
                        
                        #pkUsed = dnaAnnot.pk
                        pkUsed = lastComponentDisplayID + 1
                        stringId = '00000000'+str(pkUsed)                                
                        stringId = stringId[-3:]
                        dnaAnnot.displayId = 'Ve'+stringId
                        dnaAnnot.save()
                    else:
                        lastComponent = DnaComponent.objects.filter(displayId__regex='^gb.*').order_by('-displayId')[:1]
                        leng = len(lastComponent)
                        if leng>0:
                            lastComponentDisplayID = int(lastComponent[0].displayId[2:])
                        else:
                            lastComponentDisplayID = 0  
                        pkUsed = lastComponentDisplayID + 1    
                        stringId = '00000000'+str(pkUsed)                                
                        stringId = stringId[-4:]
                        dnaAnnot.displayId = 'gb'+stringId
                        dnaAnnot.save()                                
                dnaAnnot.componentType = subCtVectorType
                dnaAnnot.save()
                stra = utilLabrack.getStrand(fullSequence,dnaAnnot.sequence)
               
                an2db = DnaSequenceAnnotation(uri ='', bioStart = first, bioEnd = secon, strand = stra, subComponent = instance, componentAnnotated = dnaAnnot)
                an2db.save()
               
        print '+ end selectedAnnotFromGB'
        
        return instance

    class Meta:
        model = DnaComponent  