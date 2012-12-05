from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from tyers_site.labrack.models import CSVSample
from tyers_site.labrack.models import Sample
from tyers_site import settings
from django.core.files.storage import FileSystemStorage
from labrack.models.generalmodels import Container
from labrack.models.component import DnaComponent
from labrack.models.component import ChemicalComponent
from labrack.models.component import PeptideComponent
from labrack.models.component import ProteinComponent
from labrack.models.sample import SampleContent
from labrack.models.unit import Unit
from django.contrib.contenttypes.models import ContentType
from django.contrib import messages
import os

from labrack.models.generalmodels import Document
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
                        try:
                                myCSVSsample = CSVSample.import_data(data = open(gb_file2))
                        except Exception:
                                print "Oops!  That was no valid number.  Try again..." 
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
                                #Chemical Content
                                chemLine = each_line.Chemical_DisplayID.strip()
                                if (chemLine != ''):
                                        chemicalFromDB = ChemicalComponent.objects.get(displayId=chemLine)
                                        content_type = ContentType.objects.get(model='chemicalcomponent')
                                        object_id = chemicalFromDB.id
                                        solvent = each_line.Chemical_SolventBuffer
                                        concentration = each_line.Chemical_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.Chemical_Concentration_Unit)
                                        amount = each_line.Chemical_Amount
                                        amountUnit = Unit.objects.get(name=each_line.Chemical_Amount_Unit)
                                        sampleContentToDB = SampleContent(sample = sample2db,content_type = content_type,object_id = object_id, solvent=solvent,concentration = concentration,concentrationUnit = concentrationUnit, amount=amount, amountUnit=amountUnit)
                                        sampleContentToDB.save()
                                #Peptide Content
                                peptideLine = each_line.Peptide_DisplayID.strip()
                                if (peptideLine != ''):
                                        peptideFromDB = PeptideComponent.objects.get(displayId=peptideLine)
                                        content_type = ContentType.objects.get(model='peptidecomponent')
                                        object_id = peptideFromDB.id
                                        solvent = each_line.Peptide_SolventBuffer
                                        concentration = each_line.Peptide_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.Peptide_Concentration_Unit)
                                        amount = each_line.Peptide_Amount
                                        amountUnit = Unit.objects.get(name=each_line.Peptide_Amount_Unit)
                                        sampleContentToDB = SampleContent(sample = sample2db,content_type = content_type,object_id = object_id, solvent=solvent,concentration = concentration,concentrationUnit = concentrationUnit, amount=amount, amountUnit=amountUnit)
                                        sampleContentToDB.save() 
                                #Protein Content
                                proteinLine = each_line.Protein_DisplayID.strip()
                                if (proteinLine != ''):
                                        proteinFromDB = ProteinComponent.objects.get(displayId=proteinLine)
                                        content_type = ContentType.objects.get(model='proteincomponent')
                                        object_id = proteinFromDB.id
                                        solvent = each_line.Protein_SolventBuffer
                                        concentration = each_line.Protein_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.Protein_Concentration_Unit)
                                        amount = each_line.Protein_Amount
                                        amountUnit = Unit.objects.get(name=each_line.Protein_Amount_Unit)
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




def dnalist(request):
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
                        try:
                                myCSVSsample = CSVSample.import_data(data = open(gb_file2))
                        except Exception:
                                print "Oops!  That was no valid number.  Try again..." 
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
                                #Chemical Content
                                chemLine = each_line.Chemical_DisplayID.strip()
                                if (chemLine != ''):
                                        chemicalFromDB = ChemicalComponent.objects.get(displayId=chemLine)
                                        content_type = ContentType.objects.get(model='chemicalcomponent')
                                        object_id = chemicalFromDB.id
                                        solvent = each_line.Chemical_SolventBuffer
                                        concentration = each_line.Chemical_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.Chemical_Concentration_Unit)
                                        amount = each_line.Chemical_Amount
                                        amountUnit = Unit.objects.get(name=each_line.Chemical_Amount_Unit)
                                        sampleContentToDB = SampleContent(sample = sample2db,content_type = content_type,object_id = object_id, solvent=solvent,concentration = concentration,concentrationUnit = concentrationUnit, amount=amount, amountUnit=amountUnit)
                                        sampleContentToDB.save()
                                #Peptide Content
                                peptideLine = each_line.Peptide_DisplayID.strip()
                                if (peptideLine != ''):
                                        peptideFromDB = PeptideComponent.objects.get(displayId=peptideLine)
                                        content_type = ContentType.objects.get(model='peptidecomponent')
                                        object_id = peptideFromDB.id
                                        solvent = each_line.Peptide_SolventBuffer
                                        concentration = each_line.Peptide_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.Peptide_Concentration_Unit)
                                        amount = each_line.Peptide_Amount
                                        amountUnit = Unit.objects.get(name=each_line.Peptide_Amount_Unit)
                                        sampleContentToDB = SampleContent(sample = sample2db,content_type = content_type,object_id = object_id, solvent=solvent,concentration = concentration,concentrationUnit = concentrationUnit, amount=amount, amountUnit=amountUnit)
                                        sampleContentToDB.save() 
                                #Protein Content
                                proteinLine = each_line.Protein_DisplayID.strip()
                                if (proteinLine != ''):
                                        proteinFromDB = ProteinComponent.objects.get(displayId=proteinLine)
                                        content_type = ContentType.objects.get(model='proteincomponent')
                                        object_id = proteinFromDB.id
                                        solvent = each_line.Protein_SolventBuffer
                                        concentration = each_line.Protein_Concentration
                                        concentrationUnit = Unit.objects.get(name=each_line.Protein_Concentration_Unit)
                                        amount = each_line.Protein_Amount
                                        amountUnit = Unit.objects.get(name=each_line.Protein_Amount_Unit)
                                        sampleContentToDB = SampleContent(sample = sample2db,content_type = content_type,object_id = object_id, solvent=solvent,concentration = concentration,concentrationUnit = concentrationUnit, amount=amount, amountUnit=amountUnit)
                                        sampleContentToDB.save()                                        

                                    # Redirect to the document list after POST
                        #return HttpResponseRedirect(reverse('labrack.views.list'))
                        return HttpResponseRedirect('/admin/labrack/dnasample/')
        else:
                form = DocumentForm() # A empty, unbound form

        # Load documents for the list page
        documents = Document.objects.all()

        # Render list page with the documents and the form
        return render_to_response(
                'admin/labrack/dnasample/dnalist.html',
                {'documents': documents, 'form': form},
                context_instance=RequestContext(request)
        )