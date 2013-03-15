from django import forms
from django.core.exceptions import ValidationError
from labrack.models.sample import DnaSample
from labrack.models.sample import ChassisSample
from labrack.models.generalmodels import Document
from labrack.models.component import DnaComponent
from labrack.models.component import Component
from labrack.models.component import DNAComponentType
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
import os
import re

import utilLabrack



from django.forms.forms import NON_FIELD_ERRORS


class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a CSV File'
    )

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
            try:
                dna2db = DnaComponent(displayId=self.cleaned_data['DNA_Display_ID'],description=self.cleaned_data['DNA_Description'], name =self.cleaned_data['DNA_Name'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['DNA_Status'], GenBankfile = self.cleaned_data['genebank'])
                dna2db.saveSequenceWithoutAnnotations()
                #dna2db.save()
                instance.dnaConstruct = dna2db
            except:
                commit=False
        else:                    
            instance.dnaConstruct = self.cleaned_data['dnaConstruct']


        if commit:
            instance.save()
        return instance 

    class Meta:
        model = DnaSample




class ChassisSampleForm(forms.ModelForm):

    ChassisDisplayID = forms.CharField(required=False,label="Display ID")
    ChassisName = forms.CharField(required=False,label="Name")  
    #Chassis_Description = forms.CharField(widget=forms.Textarea,required=False,label="Description")
    ChassisDescription = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False,label="Description")

    description = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False)
    class Meta:
        model = ChassisSample    
        
        
        




class DnaComponentForm(forms.ModelForm):


    htmlAttribute1 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800'}),required=False,label="")
    htmlAttribute2 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800','hidden':'true'}),required=False,label="")
    htmlAttribute3 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800','hidden':'true'}),required=False,label="")
    
    
    
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(DnaComponentForm, self).__init__(*args, **kwargs)
        
    
        
    def clean(self):
        errorMessage = ""
        errors = {}        
        
        jsonCancel = self.cleaned_data['htmlAttribute2']
        dataCancel = json.loads(jsonCancel)
        try:
            gb_checked_cancel = json.loads(dataCancel["selected_annot"][2]["gb_checked_cancel"])        
            if (gb_checked_cancel=='True'):
                errorMessage = "The user has cancel the saving procedure due to annotation deletions risks!"
                errors.setdefault('',[]).append(errorMessage)                
        except:
            print ''
        
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
        compTypeName = ''        
        for s in compType:
            if s.name =='Vector Backbone':
                compTypeName = 'Vector Backbone'
        
            
        regExp = r'^[a-zA-Z][a-zA-Z]\d\d\d\d[a-zA-Z]?$'
        if (compTypeName=='Vector Backbone'):
            regExp = r'^[vV]\d\d\d?$'
            
        if not(re.match(regExp, dispId)):
            if (compTypeName=='Vector Backbone'):
                errorMessage = "ID :"+dispId+" doesn't respect the vector pattern V123, example : v001 where v or V for Vector"
                errors.setdefault('',[]).append(errorMessage)                      
            #raise forms.ValidationError("ID :"+dispId+" doesn't respect the pattern XX1234X, example : sb0001A where sb for initial and A for version")
            errorMessage = "ID :"+dispId+" doesn't respect the pattern XX1234X, example : sb0001A where sb for initial and A for version"
            errors.setdefault('',[]).append(errorMessage)                      
        
        
        cleanedData = super(DnaComponentForm, self).clean()
        
        #retrieve vector GB information
        jsonFromWeb = self.cleaned_data['htmlAttribute1']
        fullSequence = self.cleaned_data['sequence']
        data = json.loads(jsonFromWeb)
        isGbData = data["data"][1]["gb_data"]
        vectorIdFromGB_displayid = data["data"][1]["displayid"]
        vectorIdFromGB_name = data["data"][1]["name"]
        vectorIdFromGB_description = data["data"][1]["description"]
        
        #retrieve annot GB information
        jsonFromWeb2 = self.cleaned_data['htmlAttribute2']
        data2 = json.loads(jsonFromWeb2)
        selectedAnnotFromGB = json.loads(data2["selected_annot"][1]["gb_checked"])
        
        if (self.instance.id!=None):        
            r = self.instance.number_related_ParentChildAnnotations()
            if (fullSequence<>self.instance.sequence and self.instance.sequence<>"" and r>0) :
                
                print 'dooo'
        
        
        # check if the name vector exist already in the DB
        if (isGbData=='true' and  data["data"][1]["name"]!=''):
            displayid_gb=data["data"][1]["displayid"]
            regExp = r'^[vV]\d\d\d?$'
            if not(re.match(regExp, displayid_gb)):
                errorMessage = "ID :"+displayid_gb+" doesn't respect the vector pattern V123, example : v001 where v or V for Vector"
                errors.setdefault('',[]).append(errorMessage)
            
            if (DnaComponent.objects.filter(displayId=displayid_gb)):
                errorMessage = "Please choose another name for the Vector "+displayid_gb+" since it already exist in the DB!"
                errors.setdefault('',[]).append(errorMessage)
        
        # check if the named annotation exist already in the DB
        if (selectedAnnotFromGB):                
                try:                       
                    for obj in selectedAnnotFromGB:
                        
                        idAnnot=obj["text2"]
                        regExp = r'^[a-zA-Z][a-zA-Z]\d\d\d\d[a-zA-Z]?$'
                        if not(re.match(regExp, idAnnot)):
                            errorMessage = "ID Annotations Error: "+idAnnot+" doesn't respect the pattern XX1234X, example : sb0001A where sb for initial and A for version"
                            errors.setdefault('',[]).append(errorMessage)                        
                        
                        if (DnaComponent.objects.filter(displayId=idAnnot)):
                            errorMessage = "ID Annotations Error: Please reload the GB file and choose another name for the Annotation "+idAnnot+" since it already exist in the DB!"
                            errors.setdefault('',[]).append(errorMessage)
                            
                        if (dispId==idAnnot):
                            errorMessage = "ID Annotations Error: the name for the Annotation "+idAnnot+" should be different from the created DNA Part!"
                            errors.setdefault('',[]).append(errorMessage)                        
                            
                      
                except Exception, err:
                    print err
                    commit=False                                 
        if  (len(errors)>0):
            raise forms.ValidationError(errors)
        
        return cleanedData 

    def save(self, commit=True):
        instance = super(DnaComponentForm, self).save(commit=False)
        fileName = self.cleaned_data['htmlAttribute3']
        
        #settings.MEDIA_ROOT+"/"+os.path.normpath(self.GenBankfile.name)
        if (fileName!=""):
            file = open(settings.MEDIA_ROOT+"/documents/GenBank/"+fileName)
            instance.GenBankfile = File(file)
        #instance.historyDescription = self.cleaned_data['genebank'] # etc


        jsonFromWeb = self.cleaned_data['htmlAttribute1']
        
        fullSequence = self.cleaned_data['sequence']
        try:
            created_by = self.cleaned_data['created_by']
        except Exception, err:
            print err
        data = json.loads(jsonFromWeb)


        isDbData = data["data"][0]["db_data"]
        isGbData = data["data"][1]["gb_data"]



   

        jsonFromWeb2 = self.cleaned_data['htmlAttribute2']
        data2 = json.loads(jsonFromWeb2)
        selectedAnnotFromDB = json.loads(data2["selected_annot"][0]["db_checked"])
        selectedAnnotFromGB = json.loads(data2["selected_annot"][1]["gb_checked"])


        #dnaConstructVector = self.cleaned_data['dnaConstructVector']
        isVectorBackbone = True
        if selectedAnnotFromDB and selectedAnnotFromGB:
            isVectorBackbone = False

        #self.is_vector_backbone = isVectorBackbone
        #get the partytypeVector
        if (not DNAComponentType.objects.filter(name='Vector Backbone')):
            subCtType = DNAComponentType(name = 'Vector Backbone')
            subCtType.save()
        subCtVectorType = DNAComponentType.objects.filter(name='Vector Backbone')

        isGbData = data["data"][1]["gb_data"]
        vectorIdFromGB_displayid = data["data"][1]["displayid"]
        vectorIdFromGB_name = data["data"][1]["name"]
        vectorIdFromGB_description = data["data"][1]["description"]





        #if commit:
        instance.save()
        #delete all the entry for annotation for this component ID and recreate if needed
        DnaSequenceAnnotation.deleteRelated(instance.displayId)        

        ### saving Vector annotation
        if (isDbData=='true' and  data["data"][0]["name"]!='') :
            vectorIdFromDB = data["data"][0]["name"]
            dnaVector = DnaComponent.objects.get(displayId=vectorIdFromDB)
            first = fullSequence.find(dnaVector.sequence)
            
                                
            second = first + len(dnaVector.sequence)
            stra = utilLabrack.getStrand(fullSequence,dnaVector.sequence)
            an2db = DnaSequenceAnnotation(uri ='', bioStart = first, bioEnd = second, strand = stra, subComponent = instance, componentAnnotated = dnaVector)
            
            an2db.save()		    

        if (isGbData=='true' and  data["data"][1]["name"]!='') :
            coverage=data["data"][1]["coverage"]
            id=data["data"][1]["displayid"]
            coverage = coverage.split('-')
            first = int(coverage[0])
            secon = int(coverage[1])
            if (first>secon):
                fir = fullSequence[first:]
                sec = fullSequence[:secon]
                seq = fir+sec
            else:
                seq = fullSequence[first:secon]            
            name=data["data"][1]["name"]
            descrp=data["data"][1]["description"]

            dnaAnnot = DnaComponent(displayId=id,description=descrp, name =name, sequence = seq, created_by = created_by)
            dnaAnnot.save()
            dnaAnnot.componentType = subCtVectorType
            dnaAnnot.save()
            stra = utilLabrack.getStrand(fullSequence,seq)
            try:
                an2db = DnaSequenceAnnotation(uri ='', bioStart = first, bioEnd = secon, strand = stra, subComponent = instance, componentAnnotated = dnaAnnot)
                an2db.save()	
            except err:
                print err                

        try:
            #if (dna2db is None):
                #dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['description'], name =self.cleaned_data['name'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
                #dna2db.save()            
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

                except Exception, err:
                    print err
                    commit=False
            if (selectedAnnotFromGB):                
                try:                       
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
                        subCtVectorType = DNAComponentType.objects.filter(name=dnatype)	
                        dnaAnnot = DnaComponent(displayId=id,description=descrp, name =name, sequence = seq,optimizedFor = chassisOptimizedFor, created_by = created_by)
                        dnaAnnot.save()
                        dnaAnnot.componentType = subCtVectorType
                        dnaAnnot.save()
                        stra = utilLabrack.getStrand(fullSequence,dnaAnnot.sequence)
                       
                        an2db = DnaSequenceAnnotation(uri ='', bioStart = first, bioEnd = secon, strand = stra, subComponent = instance, componentAnnotated = dnaAnnot)
                        an2db.save()



                except Exception, err:
                    print err
                    commit=False                   
        except Exception, err:
            print err
            commit=False                  

        return instance

    class Meta:
        model = DnaComponent  