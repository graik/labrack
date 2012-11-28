from django import forms
from django.core.exceptions import ValidationError
from labrack.models.sample import DnaSample
from labrack.models.component import DnaComponent

from django.forms.fields import DateField, ChoiceField, MultipleChoiceField
from django.forms.widgets import RadioSelect, CheckboxSelectMultiple
from django.forms.extras.widgets import SelectDateWidget

from django.forms.forms import NON_FIELD_ERRORS


class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a CSV File'
    )

class DnaSampleForm(forms.ModelForm):
    #DNA_Construct = forms.ModelChoiceField(required=False,queryset=DnaComponent.objects.all(),help_text='Either select an existing DNA Construct or fill the dna description to create a new one')
    genebank = forms.FileField(required=False)
    DNA_Display_ID = forms.CharField(required=False,label="Display ID")
    sequence = forms.CharField(widget=forms.Textarea,required=False,label="or Sequence")
    FAVORITE_COLORS_CHOICES = (
            ("available", "available"),
            ("planning", "planning"),
            ("under construction", "under construction"),
            ("abondoned", "abondoned"),
            )
    
    DNA_Status = forms.ChoiceField(FAVORITE_COLORS_CHOICES,required=False)    
      
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        
        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(DnaSampleForm, self).__init__(*args, **kwargs)
 


    
        
    def clean(self):
        cleaned_data = super(DnaSampleForm, self).clean()
        data_DNA_Construct = cleaned_data.get("dnaConstruct")
        data_DNA_Display_ID = cleaned_data.get("DNA_Display_ID")
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
                dna2db = DnaComponent(displayId=self.cleaned_data['DNA_Display_ID'], sequence = self.cleaned_data['sequence'],status = self.cleaned_data['DNA_Status'], GenBankfile = None)
                dna2db.saveWithoutAnnotations()
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
        