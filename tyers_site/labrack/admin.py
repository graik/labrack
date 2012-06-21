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
class NucleicAcidComponentAdmin(admin.ModelAdmin):

    raw_id_fields = ('translatesTo', 'variantOf', 'componentType')  
    
    fieldsets = (
        (None,
         {'fields': (('displayId', 'name', 'uri'),
                     ('shortDescription', 'abstract'),
                     ('users', 'status'))}),
        ('Details',
         {'fields' : (('componentType', 'variantOf'),
                     ('optimizedFor', 'translatesTo'),
                     ('description',)),
          'description' : 'click magnifier to select related components from' + \
                          ' pop-up dialogs. Repeat to select more than one.'
          }),
        ('Sequence Details',
         {'fields': ('sequence', 'annotations'),
          'classes':('collapse',)}),
        )

    list_display = ('displayId', 'name', 'shortDescription',
                      'show_optimizedFor', 'show_translatesTo', 'status',
                      'isavailable_cells',
                      'isavailable_dna'
                      )
    list_filter = ('users', 'status', 'optimizedFor',
                      'componentType',
                      )
    ordering = ('displayId',)
    search_fields = ('displayId', 'name', 'shortDescription', 'description',
                      'users__username', 'sequence')




    
    

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
                         'fields' : ((('displayId', 'shortDescription'),
                                      ('container', 'aliquotnb', 'preparation_date'),
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




admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(Sample, SampleAdmin) 
admin.site.register(Unit, UnitAdmin) 



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



#TODO check if needed
admin.site.register(NucleicAcidComponent, NucleicAcidComponentAdmin)
admin.site.register(ProteinComponent)
#admin.site.register(Chassis)
admin.site.register(ComponentType)
#admin.site.register(Collection)








 



