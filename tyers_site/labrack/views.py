# Create your views here.
from labrack.models.models import DnaComponentType
from django.core import serializers
from django.http import HttpResponse
from django.template import loader, Context
from datetime import datetime
from django.template import RequestContext, loader
from labrack.models import DnaComponent 

def getTypeDnaInfo(request, maintype):
    #maintype = 'Plasmid'
    currentMainType = DnaComponentType.objects.filter(subTypeOf__name=maintype)
    
    json_models = serializers.serialize("json", currentMainType)
    return HttpResponse(json_models, mimetype="application/javascript")

def getParentTypeDnaInfo(request, subtype):
    #maintype = 'Plasmid'
    currentSubType = DnaComponentType.objects.get(id=subtype)
    currentMainType = DnaComponentType.objects.filter(id = currentSubType.subTypeOf.id)
    
    json_models = serializers.serialize("json", currentMainType)
    return HttpResponse(json_models, mimetype="application/javascript")



def reviewdna(request, displayId):
    now = datetime.now()
    dnaComponent = DnaComponent.objects.get(displayId=displayId)
    t = loader.get_template('review_form.html')
    
    html = t.render(Context({'dnaComponent': dnaComponent}))
    return HttpResponse(html)