from django import forms
from labrack.models import DnaComponent, DnaComponentType, Attachement, \
     PlasmidSample, Cell,CellType



class DnaComponentForm(forms.ModelForm):
    componentType = forms.ModelMultipleChoiceField(label='Categorie',queryset=DnaComponentType.objects.filter(subTypeOf=None),required=False)
    componentSubType = forms.ModelMultipleChoiceField(label='Type',queryset=DnaComponentType.objects.all(),required=False)
    description = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    sequence = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        super(DnaComponentForm, self).__init__(*args, **kwargs)
                
    class Meta:
        model = DnaComponent 
        
class CellForm(forms.ModelForm):
    
    componentType = forms.ModelMultipleChoiceField(label='Categorie',queryset=CellType.objects.filter(subTypeOf=None),required=False)
    componentSubType = forms.ModelMultipleChoiceField(label='Type',queryset=CellType.objects.all(),required=False)
    description = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
           
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        super(CellForm, self).__init__(*args, **kwargs)
                
    class Meta:
        model = Cell 
        
        
class PlasmidSampleForm(forms.ModelForm):
    hostGenotype = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    parentalVector = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    synonyms = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
        
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        super(PlasmidSampleForm, self).__init__(*args, **kwargs)
                
    class Meta:
        model = PlasmidSample         
        

class CellSampleForm(forms.ModelForm):
    hostGenotype = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    parentalVector = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    synonyms = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
        
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        super(PlasmidSampleForm, self).__init__(*args, **kwargs)
                
    class Meta:
        model = PlasmidSample     