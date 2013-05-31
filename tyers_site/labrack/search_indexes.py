import datetime
from haystack import indexes
from labrack.models.models import DnaComponent



class DnaComponentIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    displayId = indexes.CharField(model_attr='displayId')
    registration_date = indexes.DateTimeField(model_attr='registration_date')

    def get_model(self):
        return DnaComponent

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return self.get_model().objects.filter(registration_date__lte=datetime.datetime.now())
    
