## Copyright 2012 Raik Gruenberg / Mathieu Courcelles

## This file is part of the labhamster project (http://labhamster.sf.net)
## Labhamster is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.

## Labhamster is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public
## License along with labhamster. If not, see <http://www.gnu.org/licenses/>.


from genericcollection import GenericCollectionTabularInline
from collections import OrderedDict
from django.core.exceptions import PermissionDenied
from django.core.exceptions import ValidationError
#from django.core.urlresolvers import reverse
from django.contrib import admin
from django.contrib import messages
from django.contrib.admin.sites import site
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
from adminFilters import *
from adminWidgets import *
import importexport





class PermissionAdmin():
    # Was nessary to save owners at creation
    def save_related(self, request, form, formsets, change):
        
        
        if not self.obj.pk == None:
        
            admin.ModelAdmin.save_related(self, request, form, formsets, change)
            
            ## Save the current user as an owner
            if self.obj.writePermission(request.user) or self.obj.created_by == request.user:
                self.obj.owners.add(request.user)
                self.obj.save()
            
                ## Copy Container permission to Sample
                if isinstance(self.obj, Sample):
                    self.obj.owners = self.obj.container.owners.all()
                    self.obj.owners.add(request.user)
                    self.obj.group_read = self.obj.container.group_read.all()
                    self.obj.group_write = self.obj.container.group_write.all()
                    self.obj.save()
                
                if isinstance(self.obj, Container):
                    for sample in self.obj.samples.all():
                        sample.owners = self.obj.owners.all()
                        sample.group_read = self.obj.group_read.all()
                        sample.group_write = self.obj.group_write.all()
                        sample.save()

        
    # Save the owner of the object
    def save_model(self, request, obj, form, change):

        self.obj = obj  # store the obj for save_related method

        ## Store the obj creator
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user

            # Check that user can add sample to container
            if isinstance(obj, Sample) and not obj.container.writePermission(request.user):
                messages.error(request, '%s is not allowed to add sample to this container. Ignore message below.'  
                                   % (request.user.username))
            else:
                obj.save()
        else:
            # Check if user has right to modify
            if obj.writePermission(request.user):
                obj.save()
            else:
                    messages.error(request, '%s is not allowed to modify this record. Ignore message below.'  
                                   % (request.user.username))
                    
         

            
    # Limit view to current user based on entries permission
    # Ref: http://stackoverflow.com/questions/6310983/django-admin-specific-user-admin-content
    def queryset(self, request):
        qs = admin.ModelAdmin.queryset(self, request)
        
        # Super user can see anything        
        if request.user.is_superuser:
            return qs
        

        qsSize = 0
        import re
        match = re.search('/(?P<content_type_id>\w+)/(?P<object_id>\d+)/$', request.get_full_path())
        if hasattr(match, 'group'):
            qs = qs.filter(Q(pk=int(match.group('object_id'))))
            qsSize = qs.count()


        qs = qs.filter(Q(owners=request.user) |
                         Q(group_read=request.user.groups.all()) |
                         Q(group_write=request.user.groups.all()) 
                         ).distinct()
            

        # Raise Permission denied if filtering of permission has removed required entry
        if hasattr(match, 'group') and qsSize != qs.count():
                raise PermissionDenied

                         
        return qs

    


class ComponentAdmin(PermissionAdmin, admin.ModelAdmin):
    
    
    actions = ['make_csv']
      
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                               ('Name', 'name'),
                               ('URI','uri'),
                               ('Description','description'),
                               ('Status','status'),
                               ('Abstract','abstract'),
                               ('Created by','created_by'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ]) 
    
    fieldsets = (
                 (None, {
                         'fields': (('displayId', 'name','status',),
                                    ('uri',))
                                    
                         }
                  ),
                 ('Permission', {
                                 'classes': ('collapse',),
                                 'fields' : ((('owners'), ('group_read', 'group_write'))
                                             )
                                 }
                  ),
                 ('Details', {
                              'fields' : (('variantOf', 'abstract'),
                                          ('description',)),
          
                                          }
                  ),
                 )
    
      

    list_display = ('displayId', 'name', 'created_by', 'abstract', 'status', 'showComment', )

    list_filter = ('status', 'abstract', 'created_by')
     

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
            'fields': ('GenBankfile',),
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






class DnaComponentAdmin(ComponentAdmin):
    
    fieldsets = ComponentAdmin.fieldsets.__add__(\
        ((None, {
            'fields': ('componentType','GenBankfile',),
        }
          ),('DNA Details', 
          {'fields': (('optimizedFor', 'translatesTo', 'circular'),
                      ('sequence',),
                      ),
            'classes':('collapse',)
            }
          ),)
    )
     
    
    list_display = ComponentAdmin.list_display[:-2] + \
                   ('show_optimizedFor', 'show_translatesTo', 'number_related_samples', 'size') + \
                   ComponentAdmin.list_display[-2:] 
    
    list_filter = ComponentAdmin.list_filter.__add__(('optimizedFor',))
    
    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))

    #raw_id_fields = ComponentAdmin.raw_id_fields.__add__(('translatesTo',))
    #raw_id_fields = ('componentType', 'variantOf',)
    raw_id_fields = ('componentType', 'variantOf',)
   
   

class ContainerAdmin(PermissionAdmin, admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                               ('Name', 'name'),
                               ('Container Type','containerType'),
                               ('Rack','rack'),
                               ('Created by','created_by'),
                               ('Description','description'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])
    
    fieldsets = (
                 (None, {
                         'fields' : (('displayId', 'name'),
                                     ('containerType', 'rack'),
                                     'description',
                                     )
                         }
                  ),
                 ('Permission', {
                                 'classes': ('collapse',),
                                 'fields' : ((('owners'), ('group_read', 'group_write'))
                                             )
                                 }
                  )
                 )
    
    
    list_display = ('displayId', 'name', 'containerType', 'rack_url',
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
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Container')
  


class LocationAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                               ('Name', 'name'),
                               ('Temperature','temperature'),
                               ('Room','room'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])
    
    fields = (('displayId', 'name'), 'temperature', 'room')

    list_display = ('displayId', 'name', 'temperature', 'room')
    
    list_filter = ('room', 'temperature')
    
    save_as = True
    
    search_fields = ('displayId', 'name',)
    
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Location')

    make_csv.short_description = 'Export as CSV'
    
    

class RackAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                                   ('Name', 'name'),
                                   ('Current location','location_url'),
                                                  
                                   ])
    list_display = ('displayId', 'name', 'location_url', 
                    )    
    fields = (('displayId', 'name','current_location'))
    
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
    

    
    
    


class SampleContentInline(GenericCollectionTabularInline):
    
    extra = 1
    
    model = SampleContent
       


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
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                               ('Name', 'name'),
                               ('Location', 'container.location'),
                               ('Container', 'container'),
                               ('Number of aliquots', 'aliquotNr'),
                               ('Status', 'status'),
                               ('Reference', 'reference_status'),
                               ('Description', 'description'),
                               ('Sample content', 'strFullContent()'),
                               ('Sample pedigree', 'strProvenance()'),
                               ('Sample link','sampleLinkStr()'),
                               ('Preparation date', 'preparation_date'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])

   
    fieldsets = (
                 (None, {
                         'fields' : ((('container', 'displayId'),
                                      ('preparation_date', 'aliquotNr', 'status','reference_status'),
                                      ('description'),
                                      'sampleCollection'
                                      )
                                     )
                         }
                  ),
                 
                 )
          

    inlines = [SampleContentInline, SampleProvenanceInline]
    
    list_display   = ('showId', 'location_url', 
                      'created_by', 'preparation_date', 'showSampleType', 
                      'showMainContent', 
                      'status','reference_status', 'showComment')

    list_display_links = ('showId',)
    
    list_filter = ('created_by', ContainerListFilter, 'container__rack__current_location', 
                   'status', SampleCollectionListFilter
                   )
    
    ordering       = ('container', 'displayId')
    
    raw_id_fields = ('container', 'sampleCollection')
    
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
    
    
    def file_link(self):
        if self.file:
            return "<a href='%s'>download</a>" % (self.attachment.url,)
        else:
            return "No attachment"

    file_link.allow_tags = True
    
    
    def formfield_for_dbfield(self, db_field, **kwargs):
        
        if db_field.name in self.raw_id_fields:
            
            kwargs.pop("request", None)
            fType = db_field.rel.__class__.__name__
            if fType == "ManyToOneRel":
                kwargs['widget'] = VerboseForeignKeyRawIdWidget(db_field.rel, site)
            elif fType == "ManyToManyRel":
                kwargs['widget'] = VerboseManyToManyRawIdWidget(db_field.rel, site)
            return db_field.formfield(**kwargs)
        return super(SampleAdmin, self).formfield_for_dbfield(db_field, **kwargs)
    
    
    def location_url(self, obj):
        url = obj.container.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.location.__unicode__()))
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
        
        fields = OrderedDict( [('Unit name', 'name'),
                               ('Unit Type', 'unitType'),
                               ])
        
        return importexport.generate_csv(self, request, queryset, fields, 'Unit')

    make_csv.short_description = 'Export as CSV'





class ComponentTypeAdmin(admin.ModelAdmin):

    fieldsets = (
                 (None, {
                         'fields' : (('name','uri'),
                                      'subTypeOf',
                                     )
                         }
                  ),
                 )
    
    list_display   = ('name', 'uri')
    
    ordering = ('name',)
    
    search_fields = ('name',)





class SampleCollectionAdmin(PermissionAdmin, admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Project Name', 'name'),
                                 ('Detailled Description','description'),
                               ])
    
    
    fieldsets = (
                 (None, {
                         'fields' : (( 'name',  
                                       'description',
                                      )
                                     )
                         }
                  ),
                 ('Permission', {
                        'classes': ('collapse',),
                        'fields' : (
                                        (('owners'), ('group_read', 'group_write'))
                                    )
                                 }
                  )
                 )

    list_display = ('name', 'created_by', 'creation_date')
    
    save_as = True
    
    search_fields = ('name', 'description')
    
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Project')

    make_csv.short_description = 'Export as CSV'
    
    
    
    

class SequenceAnnotationAdmin(admin.ModelAdmin):
    
    #fieldsets = ComponentAdmin.fieldsets.__add__(\
    #    (('uri', {
    #        'fields': ('sequence',),
    #        'classes':('collapse',)
    #    }
    #      ),)
    #)
    
    search_fields = ComponentAdmin.search_fields.__add__(('uri',))





## Register admin panels
admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(Rack, RackAdmin)
admin.site.register(SampleCollection, SampleCollectionAdmin)
admin.site.register(Sample, SampleAdmin) 
admin.site.register(Unit, UnitAdmin) 
admin.site.register(DNAComponentType, ComponentTypeAdmin)
admin.site.register(ProteinComponentType, ComponentTypeAdmin)
admin.site.register(PeptideComponent, PeptideComponentAdmin)
admin.site.register(ProteinComponent, ProteinComponentAdmin)
admin.site.register(DnaComponent, DnaComponentAdmin)
admin.site.register(Component, ComponentAdmin)
admin.site.register(ChemicalComponent, ComponentAdmin)
#admin.site.register(SequenceAnnotation, SequenceAnnotationAdmin)
admin.site.register(Chassis) 
#admin.site.register(Collection)
