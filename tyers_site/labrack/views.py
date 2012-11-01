from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from tyers_site.labrack.models import CSVSample
from tyers_site.labrack.models import Sample
from tyers_site import settings
from django.core.files.storage import FileSystemStorage
from tyers_site.labrack.models import Container
from tyers_site.labrack.models import DnaComponent
from tyers_site.labrack.models import SampleContent
from tyers_site.labrack.models import Unit
from django.contrib.contenttypes.models import ContentType
from django.contrib import messages
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
                        #break                
                        #create a new sample for each line
                        for each_line in myCSVSsample:
                                containerFromDB = Container.objects.get(displayId=each_line.container)
                                sample2db = Sample(displayId=each_line.DisplayId, name = each_line.name, container = containerFromDB, status = each_line.status)
                                # aliquots number
                                if (each_line.NumberOfAliquots.strip() == ''):
                                        t = int(each_line.NumberOfAliquots.strip())
                                        sample2db.aliquotNr = int(each_line.NumberOfAliquots.strip())
                                # is the sample a reference or not
                                if ((each_line.isReference.strip() == '') or (each_line.isReference.strip() == '0')):
                                        sample2db.reference_status = False
                                else:
                                        sample2db.reference_status = True
                                # end of is the sample a reference or not
                                sample2db.save()
                                #DNA Content
                                emLine = each_line.DNA_DisplayID.strip()
                                if (emLine != ''):
                                        dnaFromDB = DnaComponent.objects.get(displayId=emLine)
                                        content_type = ContentType.objects.get(model='dnacomponent')
                                        object_id = dnaFromDB.id
                                        solvent = each_line.DNA_SolventBuffer
                                        concentration = each_line.DNA_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.DNA_Concentration_Unit)
                                        amount = each_line.DNA_Amount
                                        amountUnit = Unit.objects.get(name=each_line.DNA_Amount_Unit)
                                        sampleContentToDB = SampleContent(sample = sample2db,content_type = content_type,object_id = object_id, solvent=solvent,concentration = concentration,concentrationUnit = concentrationUnit, amount=amount, amountUnit=amountUnit)
                                        sampleContentToDB.save()

                                    # Redirect to the document list after POST
                        #return HttpResponseRedirect(reverse('labrack.views.list'))
                        return HttpResponseRedirect('/admin/labrack/sample/')
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