# forms.py
from django import forms
from labrack.models import SequenceAnnotation

class SequenceAnnotationForm(forms.ModelForm):
    class Meta:
        model = SequenceAnnotation
        exclude = ('created_by',)