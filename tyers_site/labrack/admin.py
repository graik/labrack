## Copyright 2012 Raik Gruenberg

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
import tyers_site.settings as settings


admin_root = "/admin/labrack"



################################################################################################################


class ComponentAdmin(admin.ModelAdmin):
    fieldsets = (
                 (None, {
                         'fields': (('displayId', 'shortDescription','status'),
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
                              'fields' : (('componentClassification', 'variantOf', 'abstract'),
                                          ('description',)),
          
                                          }
                  ),
                 
                 )
    raw_id_fields = ('variantOf', 'componentClassification')  
    list_display = ('displayId', 'shortDescription','status', 'created_by', 'creation_date', 'modification_date')
    list_filter = ('status', 'created_by')
    search_fields = ('displayId', 'shortDescription', 'description')
    ordering = ('displayId',)
    
    
    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
        obj.save()



class ProteinComponentAdmin(ComponentAdmin):
    fieldsets = ComponentAdmin.fieldsets.__add__((('Protein Details', {
                                                                        'fields': ('sequence', 'annotations'),
                                                                        'classes':('collapse',)
                                                                        }
                                                   ),))
    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))    


class PeptideComponentAdmin(ProteinComponentAdmin):
    fieldsets = ComponentAdmin.fieldsets.__add__((('Peptide Details', {
                                                                        'fields': ('sequence', 'annotations'),
                                                                        'classes':('collapse',)
                                                                        }
                                                   ),))

################################################################################################################
class NucleicAcidComponentAdmin(ComponentAdmin):
    fieldsets = ComponentAdmin.fieldsets.__add__((('Nucleic acid Details', {
                                                                        'fields': (('sequence', 'annotations'),
                                                                                    ('optimizedFor', 'translatesTo')
                                                                                   ),
                                                                        'classes':('collapse',)
                                                                        }
                                                   ),))    
    raw_id_fields = ComponentAdmin.raw_id_fields.__add__(('translatesTo',))
    list_display = ComponentAdmin.list_display.__add__(('show_optimizedFor', 'show_translatesTo', 'isavailable_cells','isavailable_dna'))
    list_filter = ComponentAdmin.list_filter.__add__(('optimizedFor',))
    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))


   
    

################################################################################################################

class ContainerAdmin(admin.ModelAdmin):
    ## options for list view in admin interface
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
    
    list_display = ('displayId', 'shortDescription', 'containerType', 'location_url', 'created_by', 'creation_date', 'modification_date')
    list_filter = ('containerType', 'location__displayId', 'location__room', 'location__temperature', 'created_by')
    
    ordering = ('displayId',)
    search_fields = ('displayId', 'shortDescription', 'description')
    save_as = True
    
    
    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
        obj.save()
        
    def location_url(self, obj):
        url = obj.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (admin_root, url, obj.location.__unicode__()))
    location_url.allow_tags = True
    location_url.short_description = 'Location'

        



################################################################################################################
class LocationAdmin(admin.ModelAdmin):
    list_display = ('displayId', 'shortDescription', 'temperature', 'room', 'creation_date', 'modification_date')
    fields = (('displayId', 'shortDescription'), 'temperature', 'room')
    list_filter = ('room', 'temperature')
    search_fields = ('displayId', 'shortDescription',)
    save_as = True



################################################################################################################
class SampleContentInline(GenericCollectionTabularInline):
    model = SampleContent



################################################################################################################
class SampleAdmin(admin.ModelAdmin):
    
   
    fieldsets = (
                 (None, {
                         'fields' : ((('displayId', 'shortDescription', 'preparation_date'),
                                      ('container', 'aliquotnb','empty'),
                                      'description')
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
          

    inlines = [SampleContentInline,]
    class Media:
        js = (settings.MEDIA_URL + '/js/genericcollection.js',)


    list_display   = ('displayId', 'shortDescription','location_url', 'container_url', 'created_by', 'preparation_date', 'creation_date', 'modification_date')
    list_filter    = ('container', 'container__location', 'created_by')
    list_display_links  = ('displayId',)
    ordering       = ('container', 'displayId')
    search_fields  = ('displayId', 'shortDescription', 'description', 'container__displayId', 'container__location__displayId', 'container__location__temperature', 'container__location__room')
    save_as        = True
    date_hierarchy = 'preparation_date'
    

    # Save the owner of the object
    def save_model(self, request, obj, form, change):
        if getattr(obj, 'created_by', None) is None:
            obj.created_by = request.user
        obj.save()
        
    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (admin_root, url, obj.container.__unicode__()))
    container_url.allow_tags = True
    container_url.short_description = 'Container'
    
    def location_url(self, obj):
        url = obj.container.location.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (admin_root, url, obj.container.location.__unicode__()))
    location_url.allow_tags = True
    location_url.short_description = 'Location'




################################################################################################################
class UnitAdmin(admin.ModelAdmin):
    list_display = ('name', 'type')
    list_filter  = (('type'),)
    ordering     = ('type', 'name')



################################################################################################################
class ComponentClassificationAdmin(admin.ModelAdmin):
    fieldsets = (
                 (None, {
                         'fields' : (('shortDescription','uri'),
                                      'subTypeOf',
                                     )
                         }
                  ),
                 )



################################################################################################################




admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
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








 



