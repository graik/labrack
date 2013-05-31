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


from collections import OrderedDict
from django.core.exceptions import PermissionDenied
from django.core.exceptions import ValidationError
from django.contrib import admin
from django.contrib import messages
from django.contrib.admin.sites import AdminSite, site 
from django.db.models import Q
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.utils.safestring import mark_safe
from django import forms
from django.forms.util import ErrorList
from django.shortcuts import render_to_response
from django.contrib.admin import SimpleListFilter

from django.contrib.admin.util import flatten_fieldsets

## Labrack imports
import tyers_site.settings as S
from labrack.models.generalmodels import Rack
from labrack.models.generalmodels import Container
from labrack.models.generalmodels import Location
from labrack.models.models import ExtraFile
from labrack.models.models import DnaComponent
from labrack.models.models import Source
from labrack.models.models import DnaComponentType
from labrack.models.models import PlasmidSample
from labrack.forms import DnaComponentForm




from labrack.models import utilLabrack
from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from django.contrib.sites.models import Site

 

class AutoAssignAdmin():
    
    # Save the owner of the object
    def save_model(self, request, obj, form, change):

        self.obj = obj  # store the obj for save_related method
        self.obj.created_by = request.user
        #self.obj.owners.add(request.user)
        self.obj.save()
        

class ContainerAdmin(AutoAssignAdmin, admin.ModelAdmin):

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
    
    #def get_readonly_fields(self, request, obj=None):
        #if obj is not None:
            #return self.readonly_fields + ('created_by',)
        #return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['created_by'] = request.user
        request.GET = data
        return super(ContainerAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    

class LocationAdmin(AutoAssignAdmin,admin.ModelAdmin):

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
    
    
class RackAdmin(AutoAssignAdmin,admin.ModelAdmin):

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
    
    
class ComponentAdmin(admin.ModelAdmin):


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
    
    list_display = ('displayId', 'name', 'created_by', 'showComment','status')
    list_filter = ('status', 'created_by')

    ordering = ('displayId',)

    search_fields = ('displayId', 'name', 'description')


    #def formfield_for_dbfield(self, db_field, **kwargs):

        #if db_field.name in self.raw_id_fields:

            #kwargs.pop("request", None)
            #fType = db_field.rel.__class__.__name__
            #if fType == "ManyToOneRel":
                #kwargs['widget'] = VerboseForeignKeyRawIdWidget(db_field.rel, site)
            #elif fType == "ManyToManyRel":
                #kwargs['widget'] = VerboseManyToManyRawIdWidget(db_field.rel, site)
            #return db_field.formfield(**kwargs)
        #return super(ComponentAdmin, self).formfield_for_dbfield(db_field, **kwargs)

    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Component')


    make_csv.short_description = 'Export as CSV'


class vectorBackboneFilter(SimpleListFilter):
    title = 'vector Backbone'
    parameter_name = 'vectorBackbone'

    def lookups(self, request, model_admin):
        #dnaOfBaseVector = DnaComponent.objects.filter(componentType__name='Base Vector')        
        dnaOfBaseVector = DnaComponent.objects.filter(componentSubType__subTypeOf=DnaComponentType.objects.get(name='Vector Backbone'))        
        r = [(c.id, c.displayId) for c in dnaOfBaseVector]
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset
        else:
            return queryset.filter(vectorBackbone__id__exact=self.value())
        
        
class markerFilter(SimpleListFilter):
    title = 'Marker'
    parameter_name = 'marker'

    def lookups(self, request, model_admin):
        #dnaOfMarker = DnaComponent.objects.filter(componentType__name='Marker') 
        dnaOfMarker = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Marker')        
        r = [(c.id, c.displayId) for c in dnaOfMarker]        
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset
        else:
            return queryset.filter(marker__id__exact=self.value())
      
       
class PlasmidSampleAdmin(ComponentAdmin):
    #form = DnaComponentForm
    

    readonly_fields = ('registration_date',)
    

    fieldsets = (
        (None, {
            'fields' : ((('container', 'displayId', 
                          'reference_status'),
                         ('aliquotNr',),
                         ('preparation_date','status'),
                         ('description'),
                         )
                        )
        }
         ),

      
        )

    #add short_description + 
    #DiplayId +Name +Insert (clickable) + Vector  (clickable) + Marker  (clickable) + User +Shorten Description + Status + Edit 
    
    # Fiktler by : status + vector(?) + marker
    list_display = ('diplayId_readOnly', 'name', 'insert_url','vectorBackbone_url', 'marker_url', 'created_by_user','shorten_description','status','dnaComponent_editOnly')

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))
    list_filter = ComponentAdmin.list_filter.__add__(( vectorBackboneFilter, markerFilter))
    

class DnaComponentAdmin(ComponentAdmin):
    form = DnaComponentForm
    

    readonly_fields = ('registration_date',)
    

    fieldsets = (
            (None, {
                'fields': (('displayId', 'name','status'),
                           ('componentType','componentSubType'),
                           ('circular',),
                           ('vectorBackbone','marker','insert' ),
                           ('created_by','owners'),
                           ('registration_date','source')
                            )
            }
             ),
            ('Details', {
                'fields' : (('description',),
                            ('sequence','attachements',)
                            )
                          }
             ),            
        )

    #add short_description + 
    #DiplayId +Name +Insert (clickable) + Vector  (clickable) + Marker  (clickable) + User +Shorten Description + Status + Edit 
    
    # Fiktler by : status + vector(?) + marker
    list_display = ('diplayId_readOnly', 'name', 'insert_url','vectorBackbone_url', 'marker_url', 'created_by_user','shorten_description','status','dnaComponent_editOnly')

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))
    list_filter = ComponentAdmin.list_filter.__add__(( vectorBackboneFilter, markerFilter))
    
      
        
    def shorten_description(self,obj):
        descrp = obj.description
        if (len(descrp)>30):
            descrp = descrp[0:30]+"..."
        return descrp
    shorten_description.short_description = 'description'
    shorten_description.help_text="Please use the following format:."
    
    def diplayId_readOnly(self, obj):
        return mark_safe('<a href="%s/%s/%s/">%s</a>' % (S.admin_root,'../../reviewdna', obj.displayId,obj.displayId))
        
        #return mark_safe('<a href="%s/%s/%s">%s</a>' % (S.search_root,'dnacomponent/objects', obj.id,obj.displayId))
    diplayId_readOnly.allow_tags = True    
    diplayId_readOnly.short_description = 'ID'
    
    def created_by_user(self,obj):
        return obj.created_by
    created_by_user.allow_tags = True
    created_by_user.short_description = 'User'

    def vectorBackbone_readOnly(self, obj):
        if (obj.vectorBackbone==None):
            return mark_safe('')
    vectorBackbone_readOnly.allow_tags = True
    vectorBackbone_readOnly.short_description = 'Vector'
    
    def vectorBackbone_url(self, obj):
        if (obj.vectorBackbone==None):
            return mark_safe('')
        else:
            url = obj.vectorBackbone.get_absolute_url()
            return mark_safe('<a href="../../%s">%s</a>' % (url, obj.vectorBackbone.__unicode__()))            
    vectorBackbone_url.allow_tags = True
    vectorBackbone_url.short_description = 'Vector'    
    
    def marker_url(self, obj):
        if (obj.marker==None):
            return mark_safe('')
        else:
            url = obj.marker.get_absolute_url()
            return mark_safe('<a href="../../%s">%s</a>' % (url, obj.marker.__unicode__()))
    marker_url.allow_tags = True
    marker_url.short_description = 'Marker'   
    
    def insert_url(self, obj):
        if (obj.insert==None):
            return mark_safe('')
        else:
            url = obj.insert.get_absolute_url()
            return mark_safe('<a href="../../%s">%s</a>' % (url, obj.insert.__unicode__()))
    insert_url.allow_tags = True
    insert_url.short_description = 'Insert'        

    def dnaComponent_editOnly(self, obj):
        
        return mark_safe('<a href="%s/%s/%s/?edit=1"><img src="http://icons.iconarchive.com/icons/custom-icon-design/office/16/edit-icon.png"/></a>' % (S.admin_root,'dnacomponent', obj.id))
        
    dnaComponent_editOnly.allow_tags = True    
    dnaComponent_editOnly.short_description = 'Edit'    
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'created_by':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        if db_field.name == 'vectorBackbone':
            dnaOfBaseVector = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Vector Backbone')        
           
            kwargs['queryset'] = dnaOfBaseVector
        
        if db_field.name == 'marker':
            dnaOfMarker = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Marker')        
            kwargs['queryset'] = dnaOfMarker
            
        if db_field.name == 'insert':
            #dnaOfInsert = DnaComponent.objects.filter(componentSubType__name='Insert')
            dnaOfInsert = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Insert')
            kwargs['queryset'] = dnaOfInsert
        
        return super(DnaComponentAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def formfield_for_manytomany(self, db_field, request, **kwargs):
        if db_field.name == "componentType":
            kwargs['queryset'] = DnaComponentType.objects.filter(subTypeOf=None)
        return super(DnaComponentAdmin, self).formfield_for_manytomany(
            db_field, request, **kwargs)    
    
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
    




## Register regular admin panels
admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(Rack, RackAdmin)
admin.site.register(DnaComponent, DnaComponentAdmin)
admin.site.register(DnaComponentType)
admin.site.register(ExtraFile)
admin.site.register(Source)
admin.site.register(PlasmidSample)
#admin.site.register(User)
#admin.site.register(Group)
#admin.site.register(Site)




