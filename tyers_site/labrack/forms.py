from django import forms
from django.core.exceptions import ValidationError
from labrack.models.sample import DnaSample
from labrack.models.sample import ChassisSample
from labrack.models.generalmodels import Document
from labrack.models.component import DnaComponent
from labrack.models.component import Component
from labrack.models.component import DNAComponentType
from labrack.models.generalmodels import SequenceAnnotation
from django.core.files import File

from django.forms.fields import DateField, ChoiceField, MultipleChoiceField
#from django.forms.widgets import RadioSelect, CheckboxSelectMultiple
from django.forms.extras.widgets import SelectDateWidget
from django.forms.widgets import *
import django.utils.simplejson as json
from tyers_site import settings
import os



from django.forms.forms import NON_FIELD_ERRORS


class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a CSV File'
    )

class DnaSampleForm(forms.ModelForm):
    #DNA_Construct = forms.ModelChoiceField(required=False,queryset=DnaComponent.objects.all(),help_text='Either select an existing DNA Construct or fill the dna description to create a new one')
     


    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)

        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(DnaSampleForm, self).__init__(*args, **kwargs)

    def clean(self):
        cleaned_data = super(DnaSampleForm, self).clean()
        data_DNA_Construct = cleaned_data.get("dnaConstruct")
        data_DNA_Display_ID = cleaned_data.get("DNA_Display_ID")
        data_DNA_Name = cleaned_data.get("DNA_Name")
        data_DNA_Description = cleaned_data.get("DNA_Description")
        che = cleaned_data.get("inChassis")
        if data_DNA_Display_ID=='':
        
            data_DNA_Display_ID=None
        if (data_DNA_Construct==None and data_DNA_Display_ID==None):
            raise forms.ValidationError("Either select an existing DNA Construct or fill the dna description to create a new one!")
        return cleaned_data        


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

    Chassis_Display_ID = forms.CharField(required=False,label="Display ID")
    Chassis_Name = forms.CharField(required=False,label="Name")  
    #Chassis_Description = forms.CharField(widget=forms.Textarea,required=False,label="Description")
    Chassis_Description = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False,label="Description")

    description = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False)
    class Meta:
        model = ChassisSample        


        
        
        
class DnaComponentForm(forms.ModelForm):

  

    htmlAttribute1 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800'}),required=False,label="")
    htmlAttribute2 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800','hidden':'true'}),required=False,label="")
    
    htmlAttribute3 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800','hidden':'true'}),required=False,label="")
    
 
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
        data = json.loads(jsonFromWeb)
        
        
        is_db_data = data["data"][0]["db_data"]
	is_gb_data = data["data"][1]["gb_data"]
	
		    
        
	#var index = document.getElementById('datalist_gen').selectedIndex -1;
	#alert(arrayOfJSON_AnnotationsFile[index].coverage);        
        
        
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
	if (not DNAComponentType.objects.filter(name='Vector')):
            subCtType = DNAComponentType(name = 'Vector')
            subCtType.save()  
        subCtVectorType = DNAComponentType.objects.filter(name='Vector')	

        is_gb_data = data["data"][1]["gb_data"]        
        vectorIdFromGB_name = data["data"][1]["name"]
        vectorIdFromGB_description = data["data"][1]["description"]


        #if isVectorBackbone==True:
            #try:
                #dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['description'], name =self.cleaned_data['name'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
                #dna2db.circular = True
                #dna2db.save()
                #if (not DNAComponentType.objects.filter(name='Vector')):
                    #subCtType = DNAComponentType(name = 'Vector')
                    #subCtType.save()  
                #subCtVectorType = DNAComponentType.objects.filter(name='Vector')
                #dna2db.componentType = subCtVectorType
                #dna2db.sequence = self.cleaned_data['sequence']
                #dna2db.save()                
            #except (RuntimeError, TypeError, NameError):
                #commit=False
        #else:
            # Vector saving part
              
        
            
        
        #if commit:
	instance.save()
            
        ### saving Vector annotation
	if (is_db_data=='true' and  data["data"][0]["name"]!='') :
	    vectorIdFromDB = data["data"][0]["name"]
	    dnaVector = DnaComponent.objects.get(displayId=vectorIdFromDB)
	    first = fullSequence.find(dnaVector.sequence)
	    second = first + len(dnaVector.sequence)
	     
            an2db = SequenceAnnotation(uri ='', bioStart = first, bioEnd = second, strand = '-', subComponent = instance, componentAnnotated = dnaVector)
            an2db.save()		    
		
	if (is_gb_data=='true' and  data["data"][1]["name"]!='') :
	    coverage=data["data"][1]["coverage"]
	    id=data["data"][1]["name"]
	    coverage = coverage.split('-')
	    first = int(coverage[0])
	    secon = int(coverage[1])
	    seq = fullSequence[first:secon]
	    name=data["data"][1]["name"]
	    descrp=data["data"][1]["description"]
	    
	    dnaAnnot = DnaComponent(displayId=id,description=descrp, name =name, sequence = seq)
	    dnaAnnot.save()
	    dnaAnnot.componentType = subCtVectorType
	    dnaAnnot.save()
	    an2db = SequenceAnnotation(uri ='', bioStart = first, bioEnd = secon, strand = '-', subComponent = instance, componentAnnotated = dnaAnnot)
	    an2db.save()	

        
        #if (is_gb_data=='true'):
                    #try:
                        #dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['description'], name =self.cleaned_data['name'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
                        #dna2db.circular = True                    
                        #dna2db.save()
                        #dnaVector = DnaComponent.objects.get(displayId=vectorIdFromDB)
                        #an2db = SequenceAnnotation(uri ='', bioStart = 1, bioEnd = 2, strand = '-', subComponent = dna2db, componentAnnotated = dnaVector)
                        #an2db.save()
                        
                    #except Exception, err:
                        #print err
                        #commit=False          
	try:
	    #if (dna2db is None):
		#dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['description'], name =self.cleaned_data['name'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
		#dna2db.save()
	    if (selectedAnnotFromDB):
		try:                       
		    for id in selectedAnnotFromDB:
			dnaAnnot = DnaComponent.objects.get(id=id["id"])
			coverage=id["text6"]
			coverage = coverage.split('-')
			first = int(coverage[0])
			secon = int(coverage[1])                            
			tp = instance.sequence.index(dnaAnnot.sequence)
			an2db = SequenceAnnotation(uri ='a', bioStart = first, bioEnd = secon, strand = '-', subComponent = instance, componentAnnotated = dnaAnnot)
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
			seq = fullSequence[first:secon]
			name=obj["text3"]
			descrp=obj["text4"]
			dnatype=obj["text5"]
			optimizedfor=obj["text7"]
			
			dnaAnnot = DnaComponent(displayId=id,description=descrp, name =name, sequence = seq)
			dnaAnnot.save()
			an2db = SequenceAnnPositionotation(uri ='', bioStart = first, bioEnd = secon, strand = '-', subComponent = instance, componentAnnotated = dnaAnnot)
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