from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from tyers_site.labrack.models import CSVSample
from tyers_site import settings
from django.core.files.storage import FileSystemStorage
import os

from tyers_site.labrack.models import Document
from tyers_site.labrack.forms import DocumentForm

def list(request):
    # Handle file upload
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            newdoc = Document(docfile = request.FILES['docfile'])
            newdoc.save()
            s = request.FILES['docfile']
            p = newdoc.docfile.url
            gb_file = settings.MEDIA_ROOT+os.path.normpath(newdoc.docfile.url)
            gb_file2 = settings.TEMPLATE_DIRS+os.path.normpath(newdoc.docfile.url)
            myCSVSsample = CSVSample.import_data(data = open(gb_file2))

            # Redirect to the document list after POST
            return HttpResponseRedirect(reverse('labrack.views.list'))
    else:
        form = DocumentForm() # A empty, unbound form

    # Load documents for the list page
    documents = Document.objects.all()

    # Render list page with the documents and the form
    return render_to_response(
        'admin/labrack/sample/list.html',
        {'documents': documents, 'form': form},
        context_instance=RequestContext(request)
    )
