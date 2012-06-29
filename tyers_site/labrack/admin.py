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

from tyers_site.labrack.models import *
from django.contrib import admin
from django.http import HttpResponse
from django.utils.safestring import mark_safe
from genericcollection import GenericCollectionTabularInline
from collections import OrderedDict
import tyers_site.settings as settings
import importexport
from django.db.models import Q
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse

admin_root = "/admin/labrack"



################################################################################################################


class ComponentAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Component (ID)', 'displayId'),
                               ('Short Description', 'shortDescription'),
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
                         'fields': (('displayId', 'shortDescription', 'status'),
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
                              'fields' : (('componentClassification', 'variantOf',
                                           'abstract'),
                                          ('description',)),
          
                                          }
                  ),
                 
                 )
    
    
    list_display = ('displayId', 'shortDescription', 'status', 'created_by')
    
    list_filter = ('status', 'created_by')
    
    ordering = ('displayId',)
    
    search_fields = ('displayId', 'shortDescription', 'description')
    
    raw_id_fields = ('variantOf', 'componentClassification')
    
    
    
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Component')

    make_csv.short_description = 'Export as CSV'
    
    
    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
        obj.save()


################################################################################################################
class ProteinComponentAdmin(ComponentAdmin):
    
    ComponentAdmin.exportFields['Sequence'] = 'sequence'
    
    fieldsets = ComponentAdmin.fieldsets.__add__(
                (('Protein Details', {
                                      'fields': ('sequence', 'annotations'),
                                      'classes':('collapse',)
                                      }
                  ),)
                                                 )
    
    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))    



################################################################################################################
class PeptideComponentAdmin(ProteinComponentAdmin):
    
    ComponentAdmin.exportFields['Sequence'] = 'sequence'
    
    fieldsets = ComponentAdmin.fieldsets.__add__(
                (('Peptide Details', {
                                      'fields': ('sequence', 'annotations'),
                                      'classes':('collapse',)
                                      }
                  ),)
                                                )

################################################################################################################
class NucleicAcidComponentAdmin(ComponentAdmin):
    
    ComponentAdmin.exportFields['Sequence'] = 'sequence'
    
    fieldsets = ComponentAdmin.fieldsets.__add__(
                (('Nucleic acid Details', {
                                           'fields': (('sequence', 'annotations'),
                                                      ('optimizedFor', 'translatesTo')
                                                      ),
                                           'classes':('collapse',)
                                           }
                ),)
                                                 )
        

    list_display = ComponentAdmin.list_display.__add__(('show_optimizedFor',
                                                        'show_translatesTo',
                                                        'isavailable_cells',
                                                        'isavailable_dna'))

    list_filter = ComponentAdmin.list_filter.__add__(('optimizedFor',))

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))

    raw_id_fields = ComponentAdmin.raw_id_fields.__add__(('translatesTo',))
   
    

################################################################################################################

class ContainerAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Container (ID)', 'displayId'),
                               ('Short Description', 'shortDescription'),
                               ('Container Type','containerType'),
                               ('Location','location'),
                               ('Created by','created_by'),
                               ('Description','description'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])
    
    fieldsets = (
                 (None, {
                         'fields' : (('displayId', 'shortDescription'),
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
    
    list_display = ('displayId', 'shortDescription', 'containerType',
                    'location_url', 'created_by')
    
    list_filter = ('containerType', 'location__displayId', 'location__room',
                   'location__temperature', 'created_by')
    
    ordering = ('displayId',)
    
    save_as = True
    
    search_fields = ('displayId', 'shortDescription', 'description')
    
    
    context = {}
    context['teste'] = 'fdsfsdffsd'
    
        
    def location_url(self, obj):
        url = obj.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (admin_root, url, obj.location.__unicode__()))
    
    location_url.allow_tags = True
    
    location_url.short_description = 'Location'
    
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Container')

    make_csv.short_description = 'Export as CSV'
    
    
    
    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        
        permission_notok = 1
        
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
            obj.save()
            permission_notok = 0
        
        # Superusers can modify
        if request.user.is_superuser:
            obj.save()
            permission_notok = 0
        
        # Creator can modify
        if obj.created_by == request.user:
            obj.save()
            permission_notok = 0
        
        # Owners can modify
        for user in obj.owners.all():
            if user == request.user:
                obj.save()
                permission_notok = 0
                break
            
        # Groups with write can modify
        for group_obj in obj.group_write.all():
            for group_member in request.user.groups.all():
                if group_obj == group_member:
                    obj.save()
                    permission_notok = 0
                    break
        
        if permission_notok:
            from django.contrib import messages
            messages.error(request, '%s is not allowed to modify this record. Ignore message below.'  
                          % (request.user.username))

            
        
    
    

    # Limit view to current user based on entries permission
    # Ref: http://stackoverflow.com/questions/6310983/django-admin-specific-user-admin-content
    

    

    def change_view(self, request, object_id, form_url='', extra_context=None):
        #if not self.queryset(request).filter(id=object_id).exists():
        #    return HttpResponseRedirect(reverse('admin:labrack_Containermymodel_changelist'))
        
        
        
        return super(ContainerAdmin, self).change_view(request, object_id, form_url, extra_context)

    
    def queryset(self, request):
        qs = super(ContainerAdmin, self).queryset(request)
        
        if request.user.is_superuser:
            return qs
        return qs.filter(Q(created_by=request.user) |
                         Q(owners=request.user) |
                         Q(group_read=request.user.groups.all()) |
                         Q(group_write=request.user.groups.all()) 
                         )

    
        
    
        
        

################################################################################################################
class LocationAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Location (ID)', 'displayId'),
                               ('Short Description', 'shortDescription'),
                               ('Temperature','temperature'),
                               ('Room','room'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])
    
    fields = (('displayId', 'shortDescription'), 'temperature', 'room')

    list_display = ('displayId', 'shortDescription', 'temperature', 'room')
    
    list_filter = ('room', 'temperature')
    
    save_as = True
    
    search_fields = ('displayId', 'shortDescription',)
    
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Location')

    make_csv.short_description = 'Export as CSV'
    
    




################################################################################################################
class ProjectAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    exportFields = OrderedDict( [('Project Name', 'name'),
                               ('Short Description', 'shortDescription'),
                               ('Detailled Description','description'),
                               ])
    
    
    fieldsets = (
                 (None, {
                         'fields' : ((('name', 'shortDescription', 
                                       'description'),
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

    list_display = ('name', 'shortDescription', 'created_by', 'creation_date')
    
    save_as = True
    
    search_fields = ('name', 'shortDescription', 'description')
    
    
    
    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Project')

    make_csv.short_description = 'Export as CSV'
    
    
    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
        obj.save()
    
    

################################################################################################################
class SampleLinkInline(GenericCollectionTabularInline):
    
    extra = 1
    
    model = SampleLink

################################################################################################################
class SampleContentInline(GenericCollectionTabularInline):
    
    extra = 1
    
    model = SampleContent
    
    

################################################################################################################
class SamplePedigreeInline(GenericCollectionTabularInline):
    extra = 1
    
    fk_name = 'sample'

    model = SamplePedigree

    raw_id_fields = ('sample_source',)



################################################################################################################
class SampleAdmin(admin.ModelAdmin):
   
    actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']
    
    exportFields = OrderedDict( [('Sample (ID)', 'displayId'),
                               ('Short Description', 'shortDescription'),
                               ('Location', 'container.location'),
                               ('Container', 'container'),
                               ('Number of aliquots', 'aliquotnb'),
                               ('Status', 'status'),
                               ('Description', 'description'),
                               ('Sample content', 'sampleContentStr()'),
                               ('Sample pedigree', 'samplePedigreeStr()'),
                               ('Sample link','sampleLinkStr()'),
                               ('Preparation date', 'preparation_date'),
                               ('DB Entry creation date','creation_date'),
                               ('DB Entry modification date','modification_date'),
                               ])

   
    date_hierarchy = 'preparation_date'
   
    fieldsets = (
                 (None, {
                         'fields' : ((('displayId', 'shortDescription', 
                                       'preparation_date'),
                                      ('container', 'aliquotnb', 'status', 'attachment'),
                                      'description',
                                      'project')
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
          

    inlines = [SampleContentInline, SamplePedigreeInline, SampleLinkInline]
    
    list_display = ('displayId', 'shortDescription', 'location_url', 'container_url',
                     'created_by', 'status', 'preparation_date')

    list_display_links = ('displayId',)

    list_filter = ('container', 'container__location', 'created_by', 'status', 'project')

    ordering = ('container', 'displayId')
    
    save_as = True
    
    search_fields = ('displayId', 'shortDescription', 'description', 'container__displayId', 'container__location__displayId', 'container__location__temperature', 'container__location__room')
    
    
    
    
    class Media:
        js = (settings.MEDIA_URL + '/js/genericcollection.js',)
    
    

    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (admin_root, url, obj.container.__unicode__()))
    container_url.allow_tags = True
    container_url.short_description = 'Container'
    
    
    def file_link(self):
        if self.file:
            return "<a href='%s'>download</a>" % (self.attachment.url,)
        else:
            return "No attachment"

    file_link.allow_tags = True
    
    
    def location_url(self, obj):
        url = obj.container.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (admin_root, url, obj.container.location.__unicode__()))
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


    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
        obj.save()
        
        
    def update_status(self, request, queryset, status):
        
        i = queryset.update(status=status)

        self.message_user(request, '%i samples were set to %s'  
                          % (i, status))


################################################################################################################
class UnitAdmin(admin.ModelAdmin):
    
    actions = ['make_csv']
    
    list_display = ('name', 'type')
    
    list_filter = (('type'),)
    
    ordering = ('type', 'name')



    def make_csv(self, request, queryset):
        
        fields = OrderedDict( [('Unit name', 'name'),
                               ('Type', 'type'),
                               ])
        
        return importexport.generate_csv(self, request, queryset, fields, 'Unit')

    make_csv.short_description = 'Export as CSV'


################################################################################################################
class ComponentClassificationAdmin(admin.ModelAdmin):
    
    fieldsets = (
                 (None, {
                         'fields' : (('shortDescription', 'uri'),
                                      'subTypeOf',
                                     )
                         }
                  ),
                 )



################################################################################################################




admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(Project, ProjectAdmin)
admin.site.register(Sample, SampleAdmin) 
admin.site.register(Unit, UnitAdmin) 
admin.site.register(ComponentClassification, ComponentClassificationAdmin)
admin.site.register(PeptideComponent, PeptideComponentAdmin)
admin.site.register(ProteinComponent, ProteinComponentAdmin)
admin.site.register(ChemicalComponent, ComponentAdmin)
admin.site.register(NucleicAcidComponent, NucleicAcidComponentAdmin)



################################################################################################################    
################################################################################################################    
################################################################################################################    
##
## TODO
## Raik should check what is needed from here
##
################################################################################################################    
################################################################################################################    
################################################################################################################    

#admin.site.register(Chassis)
#admin.site.register(Collection)








 



