from django import forms
from django.core.exceptions import ValidationError
from labrack.models.sample import DnaSample
from labrack.models.sample import ChassisSample
from labrack.models.sample import PlasmidSample
from labrack.models.generalmodels import Document
from labrack.models.component import DnaComponent
from labrack.models.component import DNAComponentType

from django.forms.fields import DateField, ChoiceField, MultipleChoiceField
#from django.forms.widgets import RadioSelect, CheckboxSelectMultiple
from django.forms.extras.widgets import SelectDateWidget
from django.forms.widgets import *
import django.utils.simplejson as json



from django.forms.forms import NON_FIELD_ERRORS


class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a CSV File'
    )

class DnaSampleForm(forms.ModelForm):
    #DNA_Construct = forms.ModelChoiceField(required=False,queryset=DnaComponent.objects.all(),help_text='Either select an existing DNA Construct or fill the dna description to create a new one')
    genebank = forms.FileField(required=False)
    DNA_Display_ID = forms.CharField(required=False,label="Display ID")
    DNA_Description = forms.CharField(widget=forms.Textarea,required=False,label="Description")
    DNA_Name = forms.CharField(required=False,label="Name")
    sequence = forms.CharField(widget=forms.Textarea,required=False,label="or Sequence")
    FAVORITE_COLORS_CHOICES = (
        ("available", "available"),
        ("planning", "planning"),
        ("under construction", "under construction"),
        ("abondoned", "abondoned"),
    )

    DNA_Status = forms.ChoiceField(FAVORITE_COLORS_CHOICES,required=False)    

    description = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False)




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


class PlasmidSampleForm(forms.ModelForm):

    Plasmid_Display_ID = forms.CharField(required=False,label="Display ID")
    Plasmid_Name = forms.CharField(required=False,label="Name")  
    #Chassis_Description = forms.CharField(widget=forms.Textarea,required=False,label="Description")
    Plasmid_Description = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False,label="Description")

    description = forms.CharField(widget=forms.TextInput(attrs={'size':'80'}),required=False)


    Plasmid_attribute1 = forms.CharField(widget=forms.TextInput(attrs={'size':'1800'}),required=False,label="")
    Plasmid_attribute2 = forms.CharField(widget=forms.TextInput(attrs={'size':'18800'}),required=False,label="")

    def save(self, commit=True):
        instance = super(PlasmidSampleForm, self).save(commit=False)
        #instance.historyDescription = self.cleaned_data['genebank'] # etc

        jsonFromWeb = self.cleaned_data['Plasmid_attribute1']
        data = json.loads(jsonFromWeb)
        is_db_data = data["data"][0]["db_data"]        
        vectorIdFromDB = data["data"][0]["name"]

        is_gb_data = data["data"][1]["gb_data"]        
        vectorIdFromGB_name = data["data"][1]["name"]
        vectorIdFromGB_description = data["data"][1]["description"]


        if self.cleaned_data['is_vector_backbone']==True:
            try:
                dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['Plasmid_Description'], name =self.cleaned_data['Plasmid_Name'], sequence = self.cleaned_data['sequence_text'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
                dna2db.save()
                if (not DNAComponentType.objects.filter(name='Vector')):
                    subCtType = DNAComponentType(name = 'Vector')
                    subCtType.save()  
                subCtVectorType = DNAComponentType.objects.filter(name='Vector')
                dna2db.componentType = subCtVectorType
                dna2db.sequence = self.cleaned_data['sequence_text']
                dna2db.save()                
            except (RuntimeError, TypeError, NameError):
                commit=False
        else: # in case of annotations and vector
            
            if (is_db_data=='true'):
                try:
                    dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['Plasmid_Description'], name =self.cleaned_data['Plasmid_Name'], sequence = self.cleaned_data['sequence_text'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
                    dna2db.save()
                    dnaVector = DnaComponent.objects.filter(displayId=vectorIdFromDB)
                    an2db = SequenceAnnotation(uri ='', bioStart = 1, bioEnd = 2, strand = '-', subComponent = dna2db, componentAnnotated = dnaVector)
                    an2db.save()

                    #if (not DNAComponentType.objects.filter(name='Vector')):
                    #    subCtType = DNAComponentType(name = 'Vector')
                    #    subCtType.save()  
                    #subCtVectorType = DNAComponentType.objects.filter(name='Vector')
                    #dna2db.componentType = subCtVectorType
                    #dna2db.sequence = self.cleaned_data['sequence_text']
                    #dna2db.save()                
                except (RuntimeError, TypeError, NameError):
                    commit=False
            else:
                try:
                    dna2db = DnaComponent(displayId=self.cleaned_data['displayId'],description=self.cleaned_data['Plasmid_Description'], name =self.cleaned_data['Plasmid_Name'], sequence = self.cleaned_data['sequence_text'],status = self.cleaned_data['status'], GenBankfile = self.cleaned_data['GenBankfile'])
                    dna2db.save()


                    if (not DNAComponentType.objects.filter(name='Vector')):
                        subCtType = DNAComponentType(name = 'Vector')
                        subCtType.save()
                    dnaVector2db = DnaComponent(displayId=self.cleaned_data['displayId']+'_Vector',description=self.cleaned_data['Plasmid_Description'], name =self.cleaned_data['Plasmid_Name'], sequence = self.cleaned_data['sequence_text'],status = self.cleaned_data['status'], GenBankfile = None)
                    dnaVector2db.save()                    
                    subCtVectorType = DNAComponentType.objects.filter(name='Vector')
                    dnaVector2db.componentType = subCtVectorType
                    dnaVector2db.sequence = self.cleaned_data['sequence_text']
                    dnaVector2db.save()

                    an2db = SequenceAnnotation(uri ='', bioStart = 1, bioEnd = 2, strand = '-', subComponent = dna2db, componentAnnotated = dnaVector2db)
                    an2db.save()                    
                except (RuntimeError, TypeError, NameError):
                    commit=False                



        if commit:
            instance.save()
        return instance

    class Meta:
        model = PlasmidSample   