from django import forms
from haystack.forms import SearchForm
from labrack.models import DnaComponent


class MultipleItemSearchForm(SearchForm):
    
    marker = forms.ChoiceField( choices=(),widget=forms.Select(attrs={})) 
    marker_Search = forms.ChoiceField(choices=(('include','include'),('exclude','exclude')))
    start_date = forms.DateField(required=False)
    end_date = forms.DateField(required=False)

    def search(self):
        # First, store the SearchQuerySet received from other processing.
        sqs = super(MultipleItemSearchForm, self).search()

        if not self.is_valid():
            return self.no_query_found()

        # Check to see if a start_date was chosen.
        if self.cleaned_data['start_date']:
            sqs = sqs.filter(pub_date__gte=self.cleaned_data['start_date'])

        # Check to see if an end_date was chosen.
        if self.cleaned_data['end_date']:
            sqs = sqs.filter(pub_date__lte=self.cleaned_data['end_date'])
            
        # Check to see if an end_date was chosen.
        if self.cleaned_data['marker']:
            if self.cleaned_data['marker_Search']:
                marker_Search = self.cleaned_data['marker_Search']
                if (marker_Search=='include'):
                    sqs = sqs.filter(marker__in=self.cleaned_data['marker'])            
                else:
                    sqs = sqs.filter(marker=self.cleaned_data['marker'])            
                    

        return sqs
    
    def __init__(self, *args, **kwargs):
            super(MultipleItemSearchForm, self).__init__(*args, **kwargs)
            choiceList = []
            choiceList.append(('All','All'))
            choiceListBlock = [(i.displayId,i.displayId) for i in DnaComponent.objects.filter(componentSubType__subTypeOf__name='Marker') ]
            for q in choiceListBlock:
                choiceList.append(q)
 
            self.fields['marker'].choices = choiceList 