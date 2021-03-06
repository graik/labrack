## Copyright 2012 Raik Gruenberg / Mathieu Courcelles

## This file is part of the labrack project (http://labrack.sf.net)
## labrack is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.

## labrack is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public
## License along with labrack. If not, see <http://www.gnu.org/licenses/>.


from genericcollection import GenericCollectionTabularInline
from collections import OrderedDict
from django.core.exceptions import PermissionDenied
from django.core.exceptions import ValidationError
#from django.core.urlresolvers import reverse
from django.contrib import admin
from django.contrib import messages
from django.contrib.admin.sites import AdminSite, site 
from django.db.models import Q
from django.http import HttpResponse
from django.http import HttpResponseRedirect
#from django.utils.encoding import smart_unicode
#from django.utils.html import escape
from django.utils.safestring import mark_safe
from django import forms
from django.forms.util import ErrorList
from django.shortcuts import render_to_response





#from django.utils.translation import ugettext_lazy as _


## Labrack imports
import tyers_site.settings as S
from tyers_site.labrack.models import *
from tyers_site.labrack.models.sample import *



from adminFilters import *
from adminWidgets import *
import importexport

from labrack.models.generalmodels import Rack
from labrack.models.component import ProteinComponentType
from labrack.models.component import ProteinComponent
from labrack.models.component import DnaComponent
from labrack.models.component import DnaComponentType
from labrack.models.component import ChassisComponentType
from labrack.models.component import PeptideComponent
from labrack.models.component import Component
from labrack.models.component import ChemicalComponent
from labrack.models.generalmodels import DnaSequenceAnnotation
from labrack.models.generalmodels import ProteinSequenceAnnotation

#from labrack.models.generalmodels import Collection
#from labrack.models.generalmodels import Chassis 
import utilLabrack
 


#from tyers_site.labrack.forms import DnaSampleForm
from tyers_site.labrack.forms import ChassisSampleForm
from tyers_site.labrack.forms import DnaComponentForm
from tyers_site.labrack.forms import DnaSampleForm

from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from django.contrib.sites.models import Site




class PermissionAdmin():
    
  

    #def save_related(self, request, form, formsets, change):
        #self.obj.owners.add(request.user)
        #self.obj.save()        

    # Save the owner of the object
    def save_model(self, request, obj, form, change):

        self.obj = obj  # store the obj for save_related method
        self.obj.created_by = request.user
        #self.obj.owners.add(request.user)
        self.obj.save()
        





#class ComponentAdmin(PermissionAdmin, admin.ModelAdmin):
class ComponentAdmin( admin.ModelAdmin):


    actions = ['make_csv']

 
    exportFields = OrderedDict( [('ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Description','description'),
                                 ('Status','status'),
                                 ('Created by','created_by'),
                                 ('registered at','creation_date'),
                                 ('modified at','modification_date'),
                                 ]) 

    fieldsets = (
        (None, {
            'fields': (('displayId', 'name','status',))

        }
         ),
        ('Details', {
            'fields' : (('description',)),

        }
         ),
    )
    
    list_display = ('displayId', 'name', 'show_component_type',
                    'created_by', 'showComment','status')
    list_filter = ('status', 'created_by')

    ordering = ('displayId',)

    search_fields = ('displayId', 'name', 'description')


    def formfield_for_dbfield(self, db_field, **kwargs):

        if db_field.name in self.raw_id_fields:

            kwargs.pop("request", None)
            fType = db_field.rel.__class__.__name__
            if fType == "ManyToOneRel":
                kwargs['widget'] = VerboseForeignKeyRawIdWidget(db_field.rel, site)
            elif fType == "ManyToManyRel":
                kwargs['widget'] = VerboseManyToManyRawIdWidget(db_field.rel, site)
            return db_field.formfield(**kwargs)
        return super(ComponentAdmin, self).formfield_for_dbfield(db_field, **kwargs)


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Component')


    make_csv.short_description = 'Export as CSV'



class ProteinComponentAdmin(ComponentAdmin):

    fieldsets = ComponentAdmin.fieldsets.__add__(\
        ((None, {
            'fields': ('genBankfile',),
        }
          ),
         ('Protein Details', {
             'fields': ('sequence',),
             'classes':('collapse',)
         }
          ),)
    )
    #fieldsets = ComponentAdmin.fieldsets.__add__(\
    #        (('Protein Details', {
    #            'fields': ('sequence',),
    #            'classes':('collapse',)
    #        }
    #          ),)
    #    )    

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))

    raw_id_fields = ('componentType', 'variantOf',)




class PeptideComponentAdmin(ProteinComponentAdmin):
    fieldsets = ComponentAdmin.fieldsets.__add__(\
        (('Peptide Details', {
            'fields': ('sequence',),
            'classes':('collapse',)
        }
          ),)
    )

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))
    raw_id_fields = ('componentType', 'variantOf',)


class DnaComponentTypeFilter(SimpleListFilter):
    title = 'Type'
    parameter_name = 'dnacomponent'

    def lookups(self, request, model_admin):
        dnacomponenttypes =  DnaComponentType.objects.all()
        r = [(c.id, c.name) for c in dnacomponenttypes]
        r.append((-1,'All including annotations'))
        return r

    def queryset(self, request, queryset):
        if self.value():
            val = self.value()
            if (val == '-1'):
                #return queryset.exclude(componentType__name__exact='Annotation')
                return queryset
            else:
                #DnaComponent.objects.filter(componentType__name='Vector Backbone')
                return queryset.filter(componentType__id__exact=self.value())
        else:
            return queryset.exclude(componentType__name__exact='Annotation')





class DnaComponentAdmin(ComponentAdmin):
    form = DnaComponentForm 

    readonly_fields = ('registration_date',)

    fieldsets = (
            (None, {
                'fields': (('displayId', 'name','status'),
                           ('componentType','circular', ),
                           ('created_by','owners','registration_date')
                            )
            }
             ),
            ('Details', {
                'fields' : (('variantOf','htmlAttribute1','htmlAttribute2','htmlAttribute3'),
                            ('description',),
                            ('optimizedFor', 'translatesTo',),
                            ('sequence',),
                            ('genBankfile',)
                            )
                          }
             ),            
        )


    list_display = ('displayId', 'name', 'show_component_type',
                    'show_optimizedFor',
                    'created_by', 'showComment', 'show_resistance', 
                    'size', 'number_related_samples', 'status')

    list_filter = ComponentAdmin.list_filter.__add__(('optimizedFor',DnaComponentTypeFilter))

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))

    #raw_id_fields = ComponentAdmin.raw_id_fields.__add__(('translatesTo',))
    #raw_id_fields = ('componentType', 'variantOf',)
    raw_id_fields = ('variantOf',)
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(DnaComponentAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        
        #initialUser = str(request.user)[:2]
        proposedDisplayId = utilLabrack.getNextAvailableDNAName(request.user)
        #if (proposedDisplayId==''):
        #    proposedDisplayId = initialUser+'0001a'
        data['displayId'] = proposedDisplayId

        
        request.GET = data
        return super(DnaComponentAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    
    def save_model(self, request, obj, form, change):
        # add an additional message
        if form.errorMessageForm<>"":
            messages.warning(request, form.errorMessageForm)
        super(DnaComponentAdmin, self).save_model(request, obj, form, change)    
    

class ContainerAdmin(PermissionAdmin, admin.ModelAdmin):

    actions = ['make_csv']
    
    readonly_fields = ('registration_date',)
    

    exportFields = OrderedDict( [('Container ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Container Type','containerType'),
                                 ('Rack','rack'),
                                 ('Created by','created_by'),
                                 ('Description','description'),
                                 ('Registered at','creation_date'),
                                 ('last modified','modification_date'),
                                 ])

    fieldsets = [
        (None, {
            'fields' : (('displayId', 'name'),
                        ('containerType', 'rack'),
                        'description',('created_by','owners','registration_date')                        
                        ),
          }
         ),
    ]


    list_display = ('displayId', 'name', 'containerType', 'rack_url','location_url',
                    'created_by')

    list_filter = ('containerType', 'rack__current_location__displayId', 
                   'rack__current_location__room', 
                   'rack__current_location__temperature', 'created_by')

    ordering = ('displayId',)

    search_fields = ('displayId', 'name', 'description')

    save_as = True


    def location_url(self, obj):
        url = obj.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.location.__unicode__()))

    location_url.allow_tags = True

    location_url.short_description = 'Location'

    def rack_url(self, obj):
        url = obj.rack.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.rack.__unicode__()))

    rack_url.allow_tags = True

    rack_url.short_description = 'Rack'   
    
    def location_url(self, obj):
        if obj.rack.current_location:
            url = obj.rack.current_location.get_relative_url()
        else:
            url = ''
        
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.rack.current_location))

    location_url.allow_tags = True

    location_url.short_description = 'Location'     


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Container')


    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(ContainerAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(ContainerAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    


class LocationAdmin(PermissionAdmin,admin.ModelAdmin):

    actions = ['make_csv']
    
    readonly_fields = ('registration_date',)
    

    exportFields = OrderedDict( [('Location ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Temperature','temperature'),
                                 ('Room','room'),
                                 ('Registered at','creation_date'),
                                 ('Last modified','modification_date'),
                                 ])

    fields = (('displayId', 'name'),'temperature','room','description',('created_by','owners','registration_date'))

    list_display = ('displayId', 'name', 'temperature', 'room')

    list_filter = ('room', 'temperature')

    save_as = True

    search_fields = ('displayId', 'name',)



    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Location')

    make_csv.short_description = 'Export as CSV'

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(LocationAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(LocationAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    


class RackAdmin(PermissionAdmin,admin.ModelAdmin):

    actions = ['make_csv']
    
    readonly_fields = ('registration_date',)
    

    exportFields = OrderedDict( [('Rack ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Current location','location_url'),

                                 ])
    list_display = ('displayId', 'name', 'location_url', 
                    )    
    fields = (('displayId', 'name'),'current_location','description',('created_by','owners','registration_date'))

    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.__unicode__()))

    def location_url(self, obj):
        if obj.current_location:
            url = obj.current_location.get_relative_url()
            locatAttr = obj.current_location.__unicode__()
        else:
            url = ''
            locatAttr = ''
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, locatAttr))
    location_url.allow_tags = True

    location_url.short_description = 'Location'              


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Rack')

    make_csv.short_description = 'Export as CSV'
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(RackAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(RackAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    
    





class DnaConstructInline(GenericCollectionTabularInline):
    extra = 1

    import labrack.models as M

    model = M.DnaComponent     



class SampleProvenanceInline(GenericCollectionTabularInline):
    extra = 1

    fk_name = 'sample'

    model = SampleProvenance

    raw_id_fields = ('sample_source',)




class SampleForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request', None)
        # Voila, now you can access request anywhere in your form methods by using self.request!
        super(SampleForm, self).__init__(*args, **kwargs)

    def clean(self):
        try:
            container = Container.objects.get(pk=int(self.data['container']))
        except:
            raise ValidationError('Container id is invalid')

        if not container.writePermission(self.request.user):
            self._errors['container'] = ErrorList()
            self._errors['container'].append(str(self.request.user) + 'is not allowed to add sample to this container.')
            raise ValidationError('Validation Error')

        return self.cleaned_data

    class Meta:
        model = Sample




class SampleAdmin(PermissionAdmin, admin.ModelAdmin):

    form = SampleForm     

    def get_form(self, request, obj=None, **kwargs):

        AdminForm = super(SampleAdmin, self).get_form(request, obj, **kwargs)

        class ModelFormMetaClass(AdminForm):
            def __new__(cls, *args, **kwargs):
                kwargs['request'] = request
                return AdminForm(*args, **kwargs)

        return ModelFormMetaClass   



    actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']

    date_hierarchy = 'preparation_date'

    exportFields = OrderedDict( [('Sample ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Location', 'container.location'),
                                 ('Container', 'container'),
                                 ('Number of aliquots', 'aliquotNr'),
                                 ('Status', 'status'),
                                 ('Reference', 'reference_status'),
                                 ('Description', 'description'),
                                 ('Sample content', 'strFullContent()'),
                                 ('History', 'strProvenance()'),
                                 ('Sample link','sampleLinkStr()'),
                                 ('Preparation date', 'preparation_date'),
                                 ('Last modified','modification_date'),
                                 ])


    fieldsets = (
        (None, {
            'fields' : ((('container', 'displayId', 
                          'reference_status'),
                         ('aliquotNr',),
                         ('preparation_date','status'),
                         ('description'),
                         'sampleCollection'
                         )
                        )
        }
         ),

    )


    #inlines = [SampleContentInline, SampleProvenanceInline]
    #inlines = [SampleContentInline, SampleProvenanceInline]

    list_display   = ('showId', 'location_url', 
                      'created_by', 'preparation_date',
                      'showMainContent', 
                      'status','reference_status', 'showComment')

    list_display_links = ('showId',)

    list_filter = ('created_by', ContainerListFilter, 'container__rack__current_location', 
                   'status', SampleCollectionListFilter
                   )

    ordering       = ('container', 'displayId')

    raw_id_fields = ('container',)

    save_as        = True

    save_on_top = True

    search_fields  = ('displayId', 'name', 'description', 
                      'container__displayId', 
                      'container__rack__current_location__displayId', 
                      'container__rack__current_location__location__temperature', 
                      'container__rack__current_location__location__room')



    class Media:
        js = (S.MEDIA_URL + '/js/genericcollection.js', 
              S.MEDIA_URL + '/js/list_filter_collapse.js',)

    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.__unicode__()))
    container_url.allow_tags = True
    container_url.short_description = 'Container'

    def isreference_status(self, obj):
        refStatus = ''
        if obj.reference_status:
           refStatus = 'Yes'
        return refStatus
    isreference_status.short_description = 'Reference'
    
    def file_link(self):
        if self.file:
            return "<a href='%s'>download</a>" % (self.attachment.url,)
        else:
            return "No attachment"

    file_link.allow_tags = True


     


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Sample')

    make_csv.short_description = 'Export as CSV'


    def make_bad(self, request, queryset):
        self.update_status(request, queryset, 'bad')

    make_bad.short_description = 'Mark selected entries as bad'


    def make_empty(self, request, queryset):
        self.update_status(request, queryset, 'empty')

    make_empty.short_description = 'Mark selected entries as empty'


    def make_ok(self, request, queryset):
        self.update_status(request, queryset, 'ok')


    make_ok.short_description = 'Mark selected entries as ok'

    def location_url(self, obj):
        con = obj.container
        rac = con.rack
        locat = rac.current_location
        url = locat.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.rack.current_location.__unicode__()))

    location_url.allow_tags = True
    location_url.short_description = 'Location'    

    def qr_code_img(self, obj):
        data = obj.qr_code()
        return mark_safe('<img src="http://chart.apis.google.com/chart?cht=qr&chs=55x55&chl=' + data + '" />')
        #return mark_safe(data)
    qr_code_img.allow_tags = True
    qr_code_img.short_description = 'QR code'


    def update_status(self, request, queryset, status):

        i = 0

        for obj in queryset:
            if obj.writePermission(request.user):
                obj.status = status
                obj.save()
                i += 1
            else:
                messages.error(request, '%s is not allowed to modify %s.'  
                               % (request.user.username, obj))

        self.message_user(request, '%i samples were set to %s'  
                          % (i, status))






class DnaSampleAdmin(PermissionAdmin, admin.ModelAdmin):

    form = DnaSampleForm     

    readonly_fields = ('registration_date',)
    
    def get_form(self, request, obj=None, **kwargs):

        AdminForm = super(DnaSampleAdmin, self).get_form(request, obj, **kwargs)

        class ModelFormMetaClass(AdminForm):
            def __new__(cls, *args, **kwargs):
                kwargs['request'] = request
                return AdminForm(*args, **kwargs)

        return ModelFormMetaClass   



    actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']

    date_hierarchy = 'preparation_date'

    exportFields = OrderedDict( [('Sample ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Location', 'container.location'),
                                 ('Container', 'container'),
                                 ('Number of aliquots', 'aliquotNr'),
                                 ('Status', 'status'),
                                 ('Reference', 'reference_status'),
                                 ('Description', 'description'),
                                 ('Sample content', 'strFullContent()'),
                                 ('History', 'strProvenance()'),
                                 ('Sample link','sampleLinkStr()'),
                                 ('Preparation date', 'preparation_date'),
                                 ('Registered at','creation_date'),
                                 ('Last modified','modification_date'),
                                 ])


    fieldsets = [
        (None, {
            'fields' : ((('container', 'displayId', 'status'),
                         ('preparation_date','registration_date',),
                         ('sampleCollection','reference_status'),
                         ('concentration','concentrationUnit','amount','amountUnit',),
                         ('solvent','aliquotNr',),
                         ('created_by','owners'),
                         ('description'),
                         )
                        )
        }
         ), 
        ('DNA Content',{'fields':[('dnaConstruct','inChassis')
                                  ],
                        'description': 'Either select an existing DNA Construct or fill the dna description to create a new one'                    
                        }),
        ('History',{'fields':[('derivedFrom','provenanceType','historyDescription',)],
                    'description': 'Indicate whether this sample was created from another sample.'
                    }),        
    ]


    #inlines = [DnaConstructInline]


    list_display   = ('showId', 'location_url', 
                      'created_by', 'preparation_date','content_url', 'cell_url','concentraunit', 
                      'isreference_status','status', 'showComment')

    list_display_links = ('showId',)

    list_filter = ('created_by', ContainerListFilter, 'container__rack__current_location', 
                   'status', SampleCollectionListFilter
                   )

    ordering       = ('container', 'displayId')

    raw_id_fields = ('container','dnaConstruct','inChassis','sampleCollection','derivedFrom',)

    save_as        = True

    save_on_top = True

    search_fields  = ('displayId', 'name', 'description', 
                      'container__displayId', 
                      'container__rack__current_location__displayId', 
                      'container__rack__current_location__location__temperature', 
                      'container__rack__current_location__location__room')



    class Media:
        js = (S.MEDIA_URL + '/js/genericcollection.js', 
              S.MEDIA_URL + '/js/list_filter_collapse.js',)

   
    def isreference_status(self, obj):
        refStatus = ''
        if obj.reference_status:
           refStatus = 'Yes'
        return refStatus
    isreference_status.short_description = 'Reference'  
    
    def location_url(self, obj):
        con = obj.container
        rac = con.rack
        locat = rac.current_location
        if locat:
            url = locat.get_relative_url()
            locatAttr = obj.container.rack.current_location.__unicode__()
        else:
            url = ''
            locatAttr = '' 
        #return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.rack.current_location.__unicode__()))
        return mark_safe('<a href="%s/%s">%s, %s</a>' % (S.admin_root, url, obj.container.rack.__unicode__(),locatAttr))
    location_url.allow_tags = True
    location_url.short_description = 'Location'      
    
    def content_url(self, obj):
        url = obj.dnaConstruct.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.dnaConstruct.__unicode__()))
    content_url.allow_tags = True
    content_url.short_description = 'Content'
   
        
    def cell_url(self, obj):
        urlink = ''
        if (obj.inChassis <> None):
            url = obj.inChassis.get_relative_url()
            urlink = mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.inChassis.__unicode__()))      
        return urlink
    cell_url.allow_tags = True
    cell_url.short_description = 'inCell'    
    
    def concentraunit(self,obj):
        if obj.concentration == None:
            conc = ''
        else:
            conc = str(obj.concentration)
            
        if obj.concentrationUnit == None:
            concUnit = ''
        else:
            concUnit = str(obj.concentrationUnit)        
        
        return conc + ' '+ concUnit


    

    def file_link(self):
        if self.file:
            return "<a href='%s'>download</a>" % (self.attachment.url,)
        else:
            return "No attachment"

    file_link.allow_tags = True


    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
            
        return super(DnaSampleAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        
        request.GET = data
        return super(DnaSampleAdmin, self).add_view(request, form_url="", extra_context=extra_context)

    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Sample')

    make_csv.short_description = 'Export as CSV'

    def make_bad(self, request, queryset):
        self.update_status(request, queryset, 'bad')

    make_bad.short_description = 'Mark selected entries as bad'

    def make_empty(self, request, queryset):
        self.update_status(request, queryset, 'empty')

    make_empty.short_description = 'Mark selected entries as empty'

    def make_ok(self, request, queryset):
        self.update_status(request, queryset, 'ok')


    make_ok.short_description = 'Mark selected entries as ok'



    def qr_code_img(self, obj):
        data = obj.qr_code()
        return mark_safe('<img src="http://chart.apis.google.com/chart?cht=qr&chs=55x55&chl=' + data + '" />')
        #return mark_safe(data)
    qr_code_img.allow_tags = True
    qr_code_img.short_description = 'QR code'


    def update_status(self, request, queryset, status):

        i = 0

        for obj in queryset:
            if obj.writePermission(request.user):
                obj.status = status
                obj.save()
                i += 1
            else:
                messages.error(request, '%s is not allowed to modify %s.'  
                               % (request.user.username, obj))

        self.message_user(request, '%i samples were set to %s'  
                          % (i, status))

#class ChassisAdmin(PermissionAdmin,ComponentAdmin):
class ChassisAdmin(PermissionAdmin,ComponentAdmin):
    actions = ['make_csv']
    
    readonly_fields = ('registration_date',)
    
    
    """
    exportFields = OrderedDict( [('Rack ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Current location','location_url'),
                                 ])
                                 """
    """
    formfield_overrides = {
        models.CharField: {'widget': TextInput(attrs={'size':'20'})},
        models.TextField: {'widget': Textarea(attrs={'rows':4, 'cols':40})},
    }
    """
    
    fieldsets = (
            (None, {
                'fields': (('displayId', 'name','status'),
                           ('componentType',),
                           ('created_by','owners','registration_date')
                            )
            }
             ),
            ('Details', {
                'fields' : (('variantOf',),
                            ('description',),
                            )
                          }
             ),            
        )
    

    
    list_display   = ('displayId', 'name', 'getPartType', 
                      'getVariantOf', 'created_by', 'showComment', 'status')    

    #list_filter = ComponentAdmin.list_filter.__add__(('componentType',))
     
    raw_id_fields = ('componentType', 'variantOf',)
 
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(ChassisAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(ChassisAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    
    def getPartType(self, obj):
        #vario = len(obj.variantOf.displayId
        r = obj.componentType
        s = ', '.join([a.name for a in r.all()])
        return s
    getPartType.short_description = 'Type'     
    
    def getVariantOf(self, obj):
        #vario = len(obj.variantOf.displayId
        r = obj.variantOf
        s = ', '.join([a.displayId for a in r.all()])
        return s
    getVariantOf.short_description = 'Variant Of' 

    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.__unicode__()))

    def location_url(self, obj):
        url = obj.current_location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.current_location.__unicode__()))
    location_url.allow_tags = True

    location_url.short_description = 'Location'              


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Rack')

    make_csv.short_description = 'Export as CSV'
    
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(ChassisAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(ChassisAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    
    
    
class ChassisSampleAdmin(PermissionAdmin, admin.ModelAdmin):

    form = ChassisSampleForm    
    readonly_fields = ('registration_date',)

    #actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']

    #date_hierarchy = 'preparation_date'

    


    

    fieldsets = [
        (None, {
            'fields' : ((('container', 'displayId'),
                         ('reference_status','sampleCollection'),
                         ('preparation_date','registration_date', 'status'),
                         ('solvent','concentration','concentrationUnit','amount','amountUnit','aliquotNr',),
                         ('created_by','owners'),
                         ('description')
                         
                         )
                        
                        )
        }
         ), 
        ('Cell Content',{'fields':[('chassis'),
                                  #('Chassis_Display_ID','Chassis_Name','Chassis_Description')
                                  ],
                        }),
        ('History',{'fields':[('derivedFrom','provenanceType'),
                              ('historyDescription',)]}),
    ]



    #inlines = [DnaConstructInline]


    list_display   = ('showId', 'location_url', 
                      'created_by', 'preparation_date', 'showSampleType', 
                      'showMainContent', 'isreference_status',
                      'status', 'showComment')

    list_display_links = ('showId',)

    list_filter = ('created_by', ContainerListFilter, 'container__rack__current_location', 
                   'status', SampleCollectionListFilter, 'chassis__componentType'
                   )

    ordering       = ('container', 'displayId',)

    raw_id_fields = ('container','sampleCollection','chassis','derivedFrom',)

    save_as        = True

    save_on_top = True

    search_fields  = ('displayId', 'name', 'description', 
                      'container__displayId', 
                      'container__rack__current_location__displayId', 
                      'container__rack__current_location__location__temperature', 
                      'container__rack__current_location__location__room')



    class Media:
        js = (S.MEDIA_URL + '/js/genericcollection.js', 
              S.MEDIA_URL + '/js/list_filter_collapse.js',)

    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.__unicode__()))
        
    container_url.allow_tags = True
    container_url.short_description = 'Container'
    
    def showMainContent(self, obj):
        urlink = ''
        if (obj.chassis <> None):
            url = obj.chassis.get_relative_url()
            urlink = mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.chassis.__unicode__()))      
        return urlink
    showMainContent.allow_tags = True
    showMainContent.short_description = 'Content'
    
    def isreference_status(self, obj):
        refStatus = ''
        if obj.reference_status:
           refStatus = 'Yes'
        return refStatus
    isreference_status.short_description = 'Reference'


    def file_link(self):
        if self.file:
            return "<a href='%s'>download</a>" % (self.attachment.url,)
        else:
            return "No attachment"

    file_link.allow_tags = True


    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(ChassisSampleAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def get_readonly_fields(self, request, obj=None):
        if obj is not None:
            return self.readonly_fields + ('created_by',)
        return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(ChassisSampleAdmin, self).add_view(request, form_url="", extra_context=extra_context)
   

    def location_url(self, obj):
        con = obj.container
        rac = con.rack
        locat = rac.current_location
        if locat:
            url = locat.get_relative_url()
            locatAttr = obj.container.rack.current_location.__unicode__()
        else:
            url = ''
            locatAttr = ''
        #return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.rack.current_location.__unicode__()))
        return mark_safe('<a href="%s/%s">%s, %s</a>' % (S.admin_root, url, obj.container.rack.__unicode__(),locatAttr))
    location_url.allow_tags = True
    location_url.short_description = 'Location'   


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Sample')

    make_csv.short_description = 'Export as CSV'


    def make_bad(self, request, queryset):
        self.update_status(request, queryset, 'bad')

    make_bad.short_description = 'Mark selected entries as bad'


    def make_empty(self, request, queryset):
        self.update_status(request, queryset, 'empty')

    make_empty.short_description = 'Mark selected entries as empty'


    def make_ok(self, request, queryset):
        self.update_status(request, queryset, 'ok')


    make_ok.short_description = 'Mark selected entries as ok'



    def qr_code_img(self, obj):
        data = obj.qr_code()
        return mark_safe('<img src="http://chart.apis.google.com/chart?cht=qr&chs=55x55&chl=' + data + '" />')
        #return mark_safe(data)
    qr_code_img.allow_tags = True
    qr_code_img.short_description = 'QR code'


    def update_status(self, request, queryset, status):

        i = 0

        for obj in queryset:
            if obj.writePermission(request.user):
                obj.status = status
                obj.save()
                i += 1
            else:
                messages.error(request, '%s is not allowed to modify %s.'  
                               % (request.user.username, obj))

        self.message_user(request, '%i samples were set to %s'  
                          % (i, status))



class UnitAdmin(admin.ModelAdmin):

    actions = ['make_csv']

    list_display = ('name', 'unitType')

    list_filter = (('unitType'),)

    ordering = ('unitType', 'name')



    def make_csv(self, request, queryset):

        fields = OrderedDict([('Unit name', 'name'),
                               ('Unit type', 'unitType'),
                               ])

        return importexport.generate_csv(self, request, queryset, fields, 'Unit')

    make_csv.short_description = 'Export as CSV'



class ComponentTypeAdmin(admin.ModelAdmin):

    fieldsets = (
        (None, {
            'fields' : (('name'),
                        'subTypeOf',
                        )
        }
         ),
    )

    #list_display   = ('name','subTypeOf')

    ordering = ('name',)

    search_fields = ('name',)


class SampleCollectionAdmin(PermissionAdmin, admin.ModelAdmin):

    actions = ['make_csv']

    exportFields = OrderedDict( [('Name', 'name'),
                                 ('Description','description'),
                                 ])


    fieldsets = (
        (None, {
            'fields' : (( 'name',  
                          'description',
                          )
                        )
        }
         ),
    )

    list_display = ('name', 'created_by')

    save_as = True

    search_fields = ('name', 'description')



    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Project')

    make_csv.short_description = 'Export as CSV'



class DnaSequenceAnnotationAdmin(admin.ModelAdmin):


    search_fields = ComponentAdmin.search_fields.__add__(('subComponent','componentAnnotated'))


class ProteinSequenceAnnotationAdmin(admin.ModelAdmin):


    search_fields = ComponentAdmin.search_fields.__add__(('subComponent','componentAnnotated'))



## Register regular admin panels
admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(Rack, RackAdmin)
admin.site.register(SampleCollection, SampleCollectionAdmin)
admin.site.register(Sample, SampleAdmin)
admin.site.register(Unit, UnitAdmin) 
admin.site.register(ProteinComponentType, ComponentTypeAdmin)
admin.site.register(PeptideComponent, PeptideComponentAdmin)
admin.site.register(ProteinComponent, ProteinComponentAdmin)
admin.site.register(DnaComponent, DnaComponentAdmin)
admin.site.register(DnaComponentType, ComponentTypeAdmin) 
admin.site.register(Component, ComponentAdmin)
admin.site.register(DnaSequenceAnnotation, DnaSequenceAnnotationAdmin)
admin.site.register(Chassis,ChassisAdmin) 
admin.site.register(ChassisComponentType, ComponentTypeAdmin) 
admin.site.register(Collection)
admin.site.register(ChassisSample,ChassisSampleAdmin)
admin.site.register(DnaSample,DnaSampleAdmin)
#admin.site.register(User)
#admin.site.register(Group)
#admin.site.register(Site)


