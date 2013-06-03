from django import forms
from labrack.models import DnaComponent, DnaComponentType, ExtraFile



class DnaComponentForm(forms.ModelForm):
    componentType = forms.ModelMultipleChoiceField(label='Categorie',queryset=DnaComponentType.objects.filter(subTypeOf=None),required=False)
    componentSubType = forms.ModelMultipleChoiceField(label='Type',queryset=DnaComponentType.objects.all(),required=False)
    description = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    sequence = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 4}),required=False)
    
    
    
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(DnaComponentForm, self).__init__(*args, **kwargs)
                
    class Meta:
        model = DnaComponent 