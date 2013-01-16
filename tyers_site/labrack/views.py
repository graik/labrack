from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from tyers_site.labrack.models import CSVSample
from tyers_site.labrack.models import CSVChassisSample
from tyers_site.labrack.models import CSVCell
from tyers_site.labrack.models import Sample
from tyers_site.labrack.models import Chassis
from tyers_site.labrack.models import DnaSample
from tyers_site.labrack.models import ChassisSample
from django.utils import simplejson
from tyers_site import settings
from django.core.files.storage import FileSystemStorage
from labrack.models.generalmodels import Container
from labrack.models.component import DnaComponent
from labrack.models.component import Component
from labrack.models.component import ChemicalComponent
from labrack.models.component import PeptideComponent
from labrack.models.component import ProteinComponent
from labrack.models.component import DNAComponentType
from labrack.models.sample import SampleContent
from labrack.models.unit import Unit
from django.contrib.contenttypes.models import ContentType
from django.contrib import messages
from django.shortcuts import get_object_or_404
from django.core import serializers
from django.middleware.csrf import get_token

import os
from Bio import SeqIO
from Bio.Seq import Seq


from django.views.generic.detail import BaseDetailView, \
     SingleObjectTemplateResponseMixin


import os

from labrack.models.generalmodels import Document
from tyers_site.labrack.forms import DocumentForm


from django.http import HttpResponse
from django.shortcuts import redirect
from django import forms
from django.forms.widgets import Input
import django.utils.simplejson as json



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

                                solvent = each_line.SolventBuffer  
                                concentration = each_line.Concentration
                                concentrationUnit = Unit.objects.get(name=each_line.Concentration_Unit)
                                amount = each_line.Amount
                                amountUnit = Unit.objects.get(name=each_line.Amount_Unit)                                

                                sample2db = DnaSample(amount=amount,amountUnit=amountUnit,solvent=solvent,concentration=concentration,concentrationUnit=concentrationUnit,displayId=each_line.DisplayId, container = containerFromDB, status = each_line.status)
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
                                ###########################################sample2db.save()
                                #DNA Content
                                emLine = each_line.DNA_Construct.strip()
                                if (emLine != ''):
                                        sample2db.dnaConstruct = DnaComponent.objects.get(displayId=emLine)
                                #Cell Content
                                emLine = each_line.inCell.strip()
                                if (emLine != ''):
                                        sample2db.inChassis = Chassis.objects.get(displayId=emLine)

                                try:
                                        sample2db.save()
                                        messages.add_message(request, messages.INFO, each_line.DisplayId + ' was added.')
                                except Exception,err:
                                        messages.add_message(request, messages.ERROR, each_line.DisplayId + ' was not added!!!.    ' + repr(err) )                                         



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



def celllist(request):
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
                                myCSVSsample = CSVChassisSample.import_data(data = open(gb_file2))
                        except Exception:
                                print "Oops!  That was no valid number.  Try again..." 
                        #break                
                        #create a new sample for each line
                        for each_line in myCSVSsample:
                                containerFromDB = Container.objects.get(displayId=each_line.container)

                                solvent = each_line.SolventBuffer  
                                concentration = each_line.Concentration
                                concentrationUnit = Unit.objects.get(name=each_line.Concentration_Unit)
                                amount = each_line.Amount
                                amountUnit = Unit.objects.get(name=each_line.Amount_Unit)                                

                                sample2db = ChassisSample(amount=amount,amountUnit=amountUnit,solvent=solvent,concentration=concentration,concentrationUnit=concentrationUnit,displayId=each_line.DisplayId, container = containerFromDB, status = each_line.status)
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
                                ###########################################sample2db.save()
                                #Cell Content
                                emLine = each_line.inCell.strip()
                                if (emLine != ''):
                                        sample2db.chassis = Chassis.objects.get(displayId=emLine)

                                try:
                                        sample2db.save()
                                        messages.add_message(request, messages.INFO, each_line.DisplayId + ' was added.')
                                except Exception,err:
                                        messages.add_message(request, messages.ERROR, each_line.DisplayId + ' was not added!!!.    ' + repr(err) )                                         




                                        # Redirect to the document list after POST
                        #return HttpResponseRedirect(reverse('labrack.views.list'))
                        return HttpResponseRedirect('/admin/labrack/chassissample/')
        else:
                form = DocumentForm() # A empty, unbound form

        # Load documents for the list page
        documents = Document.objects.all()

        # Render list page with the documents and the form

        return render_to_response(
                'admin/labrack/chassissample/celllist.html',
                {'documents': documents, 'form': form},
                context_instance=RequestContext(request)
        )

def cellonlylist(request):
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
                                myCSVSsample = CSVCell.import_data(data = open(gb_file2))                                
                        except Exception:
                                print "Oops!  That was no valid number.  Try again..." 
                        #break                
                        #create a new sample for each line
                        for each_line in myCSVSsample:
                                if (each_line.DisplayId != 'DisplayId'):
                                        displayId = each_line.DisplayId
                                        name = each_line.Name
                                        description	= each_line.Description                                                            


                                        sample2db = Chassis(displayId=displayId,name=name,description=description)
                                        try:
                                                sample2db.save()
                                                messages.add_message(request, messages.INFO, displayId + ' was added.')
                                        except Exception,err:
                                                messages.add_message(request, messages.ERROR, displayId + ' was not added!!!.    ' + repr(err) )                                         


                                        # Redirect to the document list after POST
                        #return HttpResponseRedirect(reverse('labrack.views.list'))

                        return HttpResponseRedirect('/admin/labrack/chassis/')
        else:
                form = DocumentForm() # A empty, unbound form

        # Load documents for the list page
        documents = Document.objects.all()

        # Render list page with the documents and the form

        return render_to_response(
                'admin/labrack/chassis/cellonlylist.html',
                {'documents': documents, 'form': form},
                context_instance=RequestContext(request)
        )




#### new test
class Html5EmailInput(Input):
        input_type = 'email'

class ContactForm(forms.Form):
        name = forms.CharField(max_length=30)
        firstname = forms.CharField(max_length=30)
        email = forms.EmailField(max_length=50, widget=Html5EmailInput())
        message = forms.CharField(max_length=1000)
        password = forms.CharField(max_length=50, widget=forms.PasswordInput())

        def clean_password(self):
                password = self.cleaned_data['password']
                length = len(password)
                if length < 8:
                        raise forms.ValidationError("Password has to be at least 8 characters long.")
                return password

def contact(request):
        if request.method == 'POST': # If the form has been submitted...
                form = ContactForm(request.POST) # A form bound to the POST data
                if form.is_valid(): # All validation rules pass
                        name = form.cleaned_data['name']
                        firstname = form.cleaned_data['firstname']
                        email = form.cleaned_data['email']
                        message = form.cleaned_data['message']
                        password = form.cleaned_data['password']
                        # do_something

                        # then return
                        if request.is_ajax():
                                return HttpResponse(content=json.dumps({'success' : '/success'}), mimetype='application/json')

                        return redirect('success') # Redirect after POST
                elif request.is_ajax():
                        errors = json.dumps(form.errors)
                        return HttpResponse(errors, mimetype='application/json')
        else:
                form = ContactForm() # An unbound form
        return render_to_response('form.html', {'form' : form}, context_instance=RequestContext(request))

def success(request):
        return HttpResponse('success')        




def plasmidlist(request):
        form = ContactForm() # An unbound form
        return render_to_response('form.html', {'form' : form}, context_instance=RequestContext(request))

def plasmidlistw(request):
        if request.is_ajax():
                if format == 'xml':
                        mimetype = 'application/xml'
                if format == 'json':
                        mimetype = 'application/javascript'
                data = serializers.serialize(format, ExampleModel.objects.all())
                #return HttpResponse(data,mimetype) 
        # If you want to prevent non XHR calls
        else:
                return HttpResponse(status=404) 
        return render_to_response(
                'admin/labrack/chassissample/celllist.html',
                {'documents': documents, 'form': form},
                context_instance=RequestContext(request)
        )    


class JSONResponseMixin(object):
        def render_to_response(self, context):
                return self.get_json_response(self.convert_context_to_json(context))
        def get_json_response(self, content, **httpresponse_kwargs):
                return HttpResponse(content, content_type='application/json', **httpresponse_kwargs)
        def convert_context_to_json(self, context):
                return simplejson.dumps(context)

class HybridDetailView(JSONResponseMixin, SingleObjectTemplateResponseMixin, BaseDetailView):
        def render_to_response(self, context):
                if self.request.is_ajax():
                        obj = context['object'].as_dict()
                        return JSONResponseMixin.render_to_response(self, obj)
                else:
                        return SingleObjectTemplateResponseMixin.render_to_response(self, context)


def search_dna_parts(request, sequence_text):   
        message = {"list_dnas": "", "extra_values": "","parttypes_values": "","optimizedfor_values": ""}
        if request.is_ajax():
                #cell = get_object_or_404(DnaComponent, sequence='MVSKGEELFTGVVPILVELDGDVN')
                #cell = get_object_or_404(DnaComponent, sequence=cell_id)
                #dnaparts = DnaComponent.objects.filter(sequence=sequence_text)
                #retrieve Vectors 100%matching
                dnapartsVectorAll = DnaComponent.objects.filter(componentType__name='Vector')
                dnaparts = [dnapart for dnapart in DnaComponent.objects.all() if (dnapart in dnapartsVectorAll and dnapart.sequence is not None and dnapart.sequence != "" and dnapart.sequence in sequence_text)]
                name = ''
                data = serializers.serialize('json', dnaparts)
                json_object = ''
                for dnapart in dnaparts :
                        id = dnapart.id
                        name = dnapart.name
                        displayid = dnapart.displayId
                        sequence = dnapart.sequence
                        description = dnapart.description
                        if json_object == '':
                                json_object = '{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+'","displayid":"'+displayid+'","description":"'+description+'"}'
                        else:
                                json_object = json_object+',{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+'","displayid":"'+displayid+'","description":"'+description+'"}'
                json_object = '['+json_object+']'
                message['list_dnas'] = json_object

                #part for retriving potential Inserts
                dnapartsInsertsAll = DnaComponent.objects.filter(componentType__name='Vector')
                dnaparts = [dnapart for dnapart in DnaComponent.objects.all() if (dnapart not in dnapartsVectorAll and dnapart.sequence is not None and dnapart.sequence != "" and dnapart.sequence in sequence_text)]
                name = ''
                data = serializers.serialize('json', dnaparts)
                json_Insertobject = ''
                for dnapart in dnaparts :
                        id = dnapart.id
                        name = dnapart.name
                        displayid = dnapart.displayId
                        sequence = dnapart.sequence
                        fisrtPeace = len(sequence)
                        secondpeace = len(sequence_text)
                        rslt = fisrtPeace/float(secondpeace)
                        coverage = str(round(rslt*100, 1) )
                        description = dnapart.description
                        if (dnapart.optimizedFor!=None):
                                optimizedFor_name = dnapart.optimizedFor.name
                        else:
                                optimizedFor_name = ''
                        if (dnapart.componentType!=None):
                                componentType_name  = ''
                                for cpType in dnapart.componentType.all():
                                        if (componentType_name==''):
                                                componentType_name = componentType_name+cpType.name
                                        else:
                                                componentType_name = componentType_name+','+cpType.name
                        else:
                                componentType_name = ''                        
                        if json_Insertobject == '':
                                json_Insertobject = '{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+'","displayid":"'+displayid+'","description":"'+description+'","coverage":"'+coverage+'","optimizedFor_name":"'+optimizedFor_name+'","componentType_name":"'+componentType_name+'"}' 
                        else:
                                json_Insertobject = json_Insertobject+',{ "id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+'","displayid":"'+displayid+'","description":"'+description+'","coverage":"'+coverage+'","optimizedFor_name":"'+optimizedFor_name+'","componentType_name":"'+componentType_name+'"}' 
                json_Insertobject = '['+json_Insertobject+']'          


                message['extra_values'] = json_Insertobject
                
                #paret retrieving all partTypes
                dnapartstypesAll = DNAComponentType.objects.all()
                json_dnaparttype = ''
                for dnaparttype in dnapartstypesAll :
                        id = dnaparttype.id
                        name = dnaparttype.name  
                        if json_dnaparttype == '':
                                json_dnaparttype = '{ "id":"'+str(id)+'","name":"'+name+'"}'
                        else:
                                json_dnaparttype = json_dnaparttype+',{ "id":"'+str(id)+'","name":"'+name+'"}'
                json_dnaparttype = '['+json_dnaparttype+']' 
                message['parttypes_values'] = json_dnaparttype
                
                #paret retrieving all partTypes
                chassisOptimizedAll = Chassis.objects.all()
                json_chassisOptimizedAll = ''
                for chas in chassisOptimizedAll :
                        id = chas.id
                        name = chas.name  
                        if json_chassisOptimizedAll == '':
                                json_chassisOptimizedAll = '{ "id":"'+str(id)+'","name":"'+name+'"}'
                        else:
                                json_chassisOptimizedAll = json_chassisOptimizedAll+',{ "id":"'+str(id)+'","name":"'+name+'"}'
                json_chassisOptimizedAll = '['+json_chassisOptimizedAll+']' 
                message['optimizedfor_values'] = json_chassisOptimizedAll                
        else:
                message = "You're the lying type, I can just tell."
        json = simplejson.dumps(message)
        return HttpResponse(json, mimetype='application/json')



def get_dna_info(request, text_value):
        message = {"list_dnas": "", "extra_values": ""}
        try:
                if request.is_ajax():
                        #cell = get_object_or_404(DnaComponent, sequence='MVSKGEELFTGVVPILVELDGDVN')
                        #cell = get_object_or_404(DnaComponent, sequence=cell_id)
                        #dnaparts = DnaComponent.objects.filter(sequence=sequence_text)
                        #dnapart = get_object_or_404(DnaComponent, displayid=text_value)
                        dnapart = DnaComponent.objects.get(displayId=text_value)
                        json_object = ''
                        id = dnapart.id
                        name = dnapart.name
                        displayid = dnapart.displayId
                        sequence = dnapart.sequence
                        description = dnapart.description
                        json_object = '[{"id":"'+str(id)+'","name":"'+name+'","sequence":"'+sequence+'","displayid":"'+displayid+'","description":"'+description+'"}]'
                        message['list_dnas'] = json_object
                        message['extra_values'] = request.method
                else:
                        message = "You're the lying type, I can just tell."
                json = simplejson.dumps(message)
                return HttpResponse(json, mimetype='application/json')
        except:
                message = 'error'
                json = simplejson.dumps(message)
                return HttpResponse(json, mimetype='application/json')                


def upload_page( request ):
        ctx = RequestContext( request, {
                'csrf_token': get_token( request ),
        } )
        return render_to_response( 'upload_page.html', ctx )


def save_upload( uploaded, filename, raw_data ):
        '''
        raw_data: if True, uploaded is an HttpRequest object with the file being
                  the raw post data
                  if False, uploaded has been submitted via the basic form
                  submission and is a regular Django UploadedFile in request.FILES
        '''
        filename = settings.MEDIA_ROOT + '/documents/GenBank/'+ filename
        try:
                from io import FileIO, BufferedWriter
                with BufferedWriter( FileIO( filename, "wb" ) ) as dest:
                        # if the "advanced" upload, read directly from the HTTP request
                        # with the Django 1.3 functionality
                        if raw_data:
                                foo = uploaded.read( 1024 )
                                while foo:
                                        dest.write( foo )
                                        foo = uploaded.read( 1024 )
                        # if not raw, it was a form upload so read in the normal Django chunks fashion
                        else:
                                for c in uploaded.chunks( ):
                                        dest.write( c )
                        # got through saving the upload, report success
                        return True
        except IOError:
                # could not open the file most likely
                pass
        return False

def ajax_upload( request ):
        if request.method == "POST":   
                if request.is_ajax( ):
                        # the file is stored raw in the request
                        upload = request
                        is_raw = True
                        # AJAX Upload will pass the filename in the querystring if it is the "advanced" ajax upload
                        try:
                                filename = request.GET[ 'qqfile' ]
                        except KeyError:
                                return HttpResponseBadRequest( "AJAX request not valid" )
                # not an ajax upload, so it was the "basic" iframe version with submission via form
                else:
                        is_raw = False
                        if len( request.FILES ) == 1:
                                # FILES is a dictionary in Django but Ajax Upload gives the uploaded file an
                                # ID based on a random number, so it cannot be guessed here in the code.
                                # Rather than editing Ajax Upload to pass the ID in the querystring,
                                # observer that each upload is a separate request,
                                # so FILES should only have one entry.
                                # Thus, we can just grab the first (and only) value in the dict.
                                upload = request.FILES.values( )[ 0 ]
                        else:
                                raise Http404( "Bad Upload" )
                        filename = upload.name

                # save the file
                success = save_upload( upload, filename, is_raw )
                genbankjson = retrieveGenBankInfo(filename)
                # let Ajax Upload know whether we saved it or not
                import json
                ret_json = { 'success': success, 'test':genbankjson,}
                return HttpResponse( json.dumps( ret_json ) )

def retrieveGenBankInfo(filename):
        gb_file = settings.MEDIA_ROOT+"/documents/GenBank/"+os.path.normpath(filename)
        gb_features = ""
        dispId = 1
        isParsingDone = False
        json_Insertobject = '';
        for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
                if (not isParsingDone):
                        for ind in xrange(len(gb_record.features)) :
                                isParsingDone = True
                                #gb_features += '\n'+ repr(gb_record.features[ind].type) + " Location start : "+ repr(gb_record.features[ind].location._start.position) + " Location end : "+ repr(gb_record.features[ind].location._end.position)
                                nameType = repr(gb_record.features[ind].type).replace("'", "")
                                strandValue = repr(gb_record.features[ind].strand)
                                startPos = repr(gb_record.features[ind].location._start.position+1)
                                endPos = repr(gb_record.features[ind].location._end.position)
                                pos = str(startPos)+'-'+str(endPos)
                                strandValue = repr(gb_record.features[ind].strand)
                                fullSequence = gb_record.seq.tostring()
                                coverage = pos
                                description = ''
                                label = repr(gb_record.features[ind].qualifiers.get('label')).replace("['","").replace("']","").replace("\\"," ")
                                if (label == 'None' ):
                                        label = repr(gb_record.features[ind].qualifiers.get('gene')).replace("['","").replace("']","").replace("\\"," ")

                                if json_Insertobject == '':
                                        json_Insertobject = '{ "id":"'+str(label)+'","name":"'+nameType+'","sequence":"'+fullSequence+'","displayid":"'+label+'","description":"'+description+'","coverage":"'+coverage+'","startPos":"'+startPos+'","endPos":"'+endPos+'","strandValue":"'+strandValue+'"}' 
                                else:
                                        json_Insertobject = json_Insertobject+',{ "id":"'+str(label)+'","name":"'+nameType+'","sequence":"'+fullSequence+'","displayid":"'+label+'","description":"'+description+'","coverage":"'+coverage+'","startPos":"'+startPos+'","endPos":"'+endPos+'","strandValue":"'+strandValue+'"}' 
        json_Insertobject = '['+json_Insertobject+']' 
        return json_Insertobject