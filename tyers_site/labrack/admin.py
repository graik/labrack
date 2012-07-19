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

from django.contrib import admin
from django.http import HttpResponse
from django.utils.safestring import mark_safe
from genericcollection import GenericCollectionTabularInline
from collections import OrderedDict
from django.db.models import Q
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages
from django.utils.translation import ugettext_lazy as _
from django.contrib.admin import SimpleListFilter
from django import forms
from django.contrib.admin.sites import site
from django.contrib.admin.widgets import ManyToManyRawIdWidget, ForeignKeyRawIdWidget
from django.utils.encoding import smart_unicode
from django.utils.html import escape
from tyers_site.labrack.models import *
import importexport
import tyers_site.settings as S





class PermissionAdmin():
    # Was nessary to save owners at creation
    def save_related(self, request, form, formsets, change):
        admin.ModelAdmin.save_related(self, request, form, formsets, change)
        
        ## Save the current user as an owner if there is no owner
        if self.obj.owners.count() == 0:
            self.obj.owners.add(request.user)
            self.obj.save()

        
    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        
        self.obj = obj  # store the obj for save_related method

        ## Store the obj creator
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
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
        
        if request.user.is_superuser:
            return qs
        return qs.filter(Q(owners=request.user) |
                         Q(group_read=request.user.groups.all()) |
                         Q(group_write=request.user.groups.all()) 
                         ).distinct()




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
                         'fields': (('displayId', 'name','status'),
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
                              'fields' : (('componentType', 'variantOf', 'abstract'),
                                          ('description',)),
          
                                          }
                  ),
                 )
    
      

    list_display = ('displayId', 'name','status', 'created_by')

    list_filter = ('status', 'created_by')

    ordering = ('displayId',)
    
    search_fields = ('displayId', 'name', 'description')
    
    raw_id_fields = ('variantOf', 'componentType')
        
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Component')


    make_csv.short_description = 'Export as CSV'






class ProteinComponentAdmin(ComponentAdmin):
    
    fieldsets = ComponentAdmin.fieldsets.__add__(\
        (('Protein Details', {
            'fields': ('sequence', 'annotations'),
            'classes':('collapse',)
        }
          ),)
    )
    
    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))    





class PeptideComponentAdmin(ProteinComponentAdmin):
    pass





class DnaComponentAdmin(ComponentAdmin):
    
    fieldsets = ComponentAdmin.fieldsets.__add__(\
        (('DNA Details', {
            'fields': (
                ('optimizedFor', 'translatesTo'),
                ('sequence', 'annotations'),
                       ),
            'classes':('collapse',)
        }
          ),)
    )
    
    list_display = ComponentAdmin.list_display.__add__(('show_optimizedFor', 
                                                        'show_translatesTo', 
##                                                        'isavailable_cells',
##                                                        'isavailable_dna',
                                                        ))
    
    list_filter = ComponentAdmin.list_filter.__add__(('optimizedFor',))
    
    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))

    raw_id_fields = ComponentAdmin.raw_id_fields.__add__(('translatesTo',))
   



class ContainerAdmin(PermissionAdmin, admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                               ('Name', 'name'),
                               ('Container Type','containerType'),
                               ('Location','location'),
                               ('Created by','created_by'),
                               ('Description','description'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])
    
    fieldsets = (
                 (None, {
                         'fields' : (('displayId', 'name'),
                                     ('containerType', 'location'),
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
    
    list_display = ('displayId', 'name', 'containerType', 'location_url', 
                    'created_by')
    
    list_filter = ('containerType', 'location__displayId', 'location__room', 
                   'location__temperature', 'created_by')
    
    ordering = ('displayId',)
    
    search_fields = ('displayId', 'name', 'description')
 
    save_as = True
    
        
    def location_url(self, obj):
        url = obj.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.location.__unicode__()))
    
    location_url.allow_tags = True
    
    location_url.short_description = 'Location'
    
    
    
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
    
    



class SampleLinkInline(GenericCollectionTabularInline):
    
    extra = 1
    
    model = SampleLink





class SampleContentInline(GenericCollectionTabularInline):
    
    extra = 1
    
    model = SampleContent
       


class SampleProvenanceInline(GenericCollectionTabularInline):
    extra = 1
    
    fk_name = 'sample'

    model = SampleProvenance

    raw_id_fields = ('sample_source',)






class ContainerListFilter(SimpleListFilter):
    # Human-readable title which will be displayed in the
    # right admin sidebar just above the filter options.
    title = _('container')

    # Parameter for the filter that will be used in the URL query.
    parameter_name = 'container__id__exact'

    def lookups(self, request, model_admin):
        """
        Returns a list of tuples. The first element in each
        tuple is the coded value for the option that will
        appear in the URL query. The second element is the
        human-readable name for the option that will appear
        in the right sidebar.
        """
        
        qs = model_admin.queryset(request)
        
        containerList = []
        containerDict = {} 
        
        for sample in qs:
            if containerDict.has_key(sample.container.pk) != True:
                containerList.append(( str(sample.container.pk), _(str(sample.container)) ))
            containerDict[sample.container.pk] = True

        return containerList


    def queryset(self, request, qs):
        
        if request.user.is_superuser:
            pass
        else:
            qs = qs.filter(Q(owners=request.user) |
                           Q(group_read=request.user.groups.all()) |
                           Q(group_write=request.user.groups.all()) 
                           ).distinct()
        
        if(self.value() == None):
            return qs
        else:
            return qs.filter(container__id__exact=self.value())
        



class SampleAdmin(PermissionAdmin, admin.ModelAdmin):
   
    actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']
    
    date_hierarchy = 'preparation_date'
    
    exportFields = OrderedDict( [('Display ID', 'displayId'),
                               ('Name', 'name'),
                               ('Location', 'container.location'),
                               ('Container', 'container'),
                               ('Number of aliquots', 'aliquotNr'),
                               ('Status', 'status'),
                               ('Description', 'description'),
                               ('Sample content', 'sampleContentStr()'),
                               ('Sample pedigree', 'samplePedigreeStr()'),
                               ('Sample link','sampleLinkStr()'),
                               ('Preparation date', 'preparation_date'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])

   
    fieldsets = (
                 (None, {
                         'fields' : ((('container', 'displayId'),
                                      ('preparation_date', 'aliquotNr', 'status'),
                                      ('description','attachment'),
                                      'sampleCollection'
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
          

    inlines = [SampleContentInline, SampleProvenanceInline, SampleLinkInline]
    
    list_display   = ('showId', 'location_url', 
                      'created_by', 'preparation_date', 'showSampleType', 
                      'showMainContent', 
                      'status', 'showComment')

    list_display_links = ('showId',)
    
    list_filter = ('created_by', ContainerListFilter, 'container__location', 
                   'status', 'sampleCollection'
                   )
    
    ordering       = ('container', 'displayId')
    
    raw_id_fields = ('container', 'sampleCollection')
    
    save_as        = True
    
    search_fields  = ('displayId', 'name', 'description', 
                      'container__displayId', 'container__location__displayId', 
                      'container__location__temperature', 
                      'container__location__room')
    
    
    
    

    class Media:
        js = (S.MEDIA_URL + '/js/genericcollection.js',)

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
        
        i = queryset.update(status=status)

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
    
    
    
    



class VerboseForeignKeyRawIdWidget(ForeignKeyRawIdWidget):
    def label_for_value(self, value):
        key = self.rel.get_related_field().name
        try:
            obj = self.rel.to._default_manager.using(self.db).get(**{key: value})
            change_url = reverse(
                "admin:%s_%s_change" % (obj._meta.app_label, obj._meta.object_name.lower()),
                args=(obj.pk,)
            )
            return '&nbsp;<strong><a href="%s">%s</a></strong>' % (change_url, escape(obj))
        except (ValueError, self.rel.to.DoesNotExist):
            return '???'

class VerboseManyToManyRawIdWidget(ManyToManyRawIdWidget):
    def label_for_value(self, value):
        values = value.split(',')
        str_values = []
        key = self.rel.get_related_field().name
        for v in values:
            try:
                obj = self.rel.to._default_manager.using(self.db).get(**{key: v})
                x = smart_unicode(obj)
                change_url = reverse(
                    "admin:%s_%s_change" % (obj._meta.app_label, obj._meta.object_name.lower()),
                    args=(obj.pk,)
                )
                str_values += ['<strong><a href="%s">%s</a></strong>' % (change_url, escape(x))]
            except self.rel.to.DoesNotExist:
                str_values += [u'???']
        return u', '.join(str_values)




    


admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(SampleCollection, SampleCollectionAdmin)
admin.site.register(Sample, SampleAdmin) 
admin.site.register(Unit, UnitAdmin) 
admin.site.register(ComponentType, ComponentTypeAdmin)
admin.site.register(PeptideComponent, PeptideComponentAdmin)
admin.site.register(ProteinComponent, ProteinComponentAdmin)
admin.site.register(DnaComponent, DnaComponentAdmin)
admin.site.register(ChemicalComponent, ComponentAdmin)
admin.site.register(Chassis)

#admin.site.register(Collection)








 



