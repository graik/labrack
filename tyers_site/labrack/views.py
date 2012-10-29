# views.py
from django.views.generic.edit import CreateView, UpdateView, DeleteView
from django.core.urlresolvers import reverse_lazy
from labrack.models import SequenceAnnotation
from labrack.forms import SequenceAnnotationForm
from django.http import HttpResponse

def index(request):
    return HttpResponse("Hello, world. You're at the poll index.")

class SequenceAnnotationCreate(CreateView):
    #model = SequenceAnnotation
    form_class = SequenceAnnotationForm
    model = SequenceAnnotation
   
    def form_valid(self, form):
        form.instance.created_by = self.request.user
        return super(SequenceAnnotationCreate, self).form_valid(form)    

class SequenceAnnotationUpdate(UpdateView):
    model = SequenceAnnotation

class SequenceAnnotationDelete(DeleteView):
    model = SequenceAnnotation
    success_url = reverse_lazy('sequenceannotation-list')