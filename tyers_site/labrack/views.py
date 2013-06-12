# Create your views here.
from labrack.models.models import DnaComponentType,CellType
from django.core import serializers
from django.http import HttpResponse
from django.template import loader, Context
from datetime import datetime
from django.template import RequestContext, loader
from labrack.models import DnaComponent ,PlasmidSample


#### set of methods used to populate categorie and type under Dna and Cell selection,
#### called by Javascript when the user select a categorie

def getTypeDnaInfo(request, maintype):
    currentMainType = DnaComponentType.objects.filter(subTypeOf__name=maintype)
    
    json_models = serializers.serialize("json", currentMainType)
    return HttpResponse(json_models, mimetype="application/javascript") 

def getTypeCellInfo(request, maintype):
    currentMainType = CellType.objects.filter(subTypeOf__name=maintype)
    
    json_models = serializers.serialize("json", currentMainType)
    return HttpResponse(json_models, mimetype="application/javascript") 

def getParentTypeDnaInfo(request, subtype):
    currentSubType = DnaComponentType.objects.get(id=subtype)
    currentMainType = DnaComponentType.objects.filter(id = currentSubType.subTypeOf.id)
    
    json_models = serializers.serialize("json", currentMainType)
    return HttpResponse(json_models, mimetype="application/javascript")  

def getParentTypeCellInfo(request, subtype):
    currentSubType = CellType.objects.get(id=subtype)
    currentMainType = CellType.objects.filter(id = currentSubType.subTypeOf.id)
    
    json_models = serializers.serialize("json", currentMainType)
    return HttpResponse(json_models, mimetype="application/javascript") 

####


### called by DNA Component View template
def reviewdna(request, displayId):
    now = datetime.now()
    dnaComponent = DnaComponent.objects.get(displayId=displayId)
    t = loader.get_template('admin/labrack/dnacomponent/review_form.html')
    
    html = t.render(Context({'dnaComponent': dnaComponent}))
    return HttpResponse(html)

### called by Plasmid View Template
def reviewplasmid(request, id):
    now = datetime.now()
    plasmidSample = PlasmidSample.objects.get(id=id)
    t = loader.get_template('admin/labrack/plasmidsample/review_form.html')
    
    html = t.render(Context({'plasmidSample': plasmidSample}))
    return HttpResponse(html)