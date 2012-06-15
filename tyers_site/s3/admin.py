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

from tyers_site.s3.models import *
from django.contrib import admin
from django.http import HttpResponse

class DnaComponentAdmin(admin.ModelAdmin):

    raw_id_fields = ('translatesTo','variantOf','componentType')  
    
    fieldsets = (
        (None,
         {'fields': (('displayId', 'name', 'uri'),
                     ('shortDescription','abstract'),
                     ('users','status'))}),
        ('Details',
         {'fields' : (('componentType','variantOf'),
                     ('optimizedFor', 'translatesTo'),
                     ('description',)),
          'description' : 'click magnifier to select related components from'+\
                          ' pop-up dialogs. Repeat to select more than one.'
          }),
        ('Sequence Details',
         {'fields': ('sequence','annotations'),
          'classes':('collapse',)} ),
        )

    list_display   = ('displayId','name','shortDescription',
                      'show_optimizedFor','show_translatesTo', 'status',
                      'isavailable_cells',
                      'isavailable_dna'
                      )
    list_filter    = ('users', 'status', 'optimizedFor',
                      'componentType',
                      )
    ordering       = ('displayId',)
    search_fields  = ('displayId','name','shortDescription','description',
                      'users__username','sequence' )

admin.site.register(DnaComponent, DnaComponentAdmin)


admin.site.register(ProteinComponent)
admin.site.register(Chassis)
admin.site.register(ComponentType)
admin.site.register(Location)


class SampleAdmin(admin.ModelAdmin):
    
    raw_id_fields = ('dna','cell','vector','protein') 
    readonly_fields = ('created',)
    
    fieldsets = (
        (None, {
            'fields' : (
                (),
                ('displayId','container','vesselType','created')
            )}),
        ('Content', {
            'fields' : (
                ('dna', 'vector', 'cell', 'protein'),
                ('concentration', 'concentrationUnit')
                ),
            'description' : 'Select sample content by clicking on magnifiers. '
            }),
         ('Additional details', {
             'fields' : (
                 ('users','comments',)
             )})
          )

    list_display   = ('show_Id', 'show_sampleType', 
                      'show_dna','show_vector','show_cell','show_protein',
                      'show_concentration', 'created', 'show_comments' )
    list_filter    = ('vesselType','container','users')
    ordering       = ('container','displayId',)
    search_fields  = ('displayId','comments','users__username',
                      'container__displayId', 
                      'dna__displayId', 'vector__displayId', 'cell__displayId',
                      'protein__displayId',
                      'dna__name', 'vector__name', 'cell__name', 
                      'protein__name')
    
    save_as        = True
    
admin.site.register(Sample, SampleAdmin)


class SampleContainerAdmin(admin.ModelAdmin):
    ## options for list view in admin interface
    fieldsets = (
        (None, {
        'fields' : (('displayId', 'shortDescription'), 
                    ('containerType','users'),
                    'description')
        }),
    )
    
    list_display   = ('displayId', 'containerType',
                      'shortDescription')
    
    list_filter    = ('containerType','users')
    ordering       = ('displayId',)
    search_fields  = ('displayId','shortDescription','description',
                      'users__username')
    save_as        = True
    
admin.site.register(SampleContainer, SampleContainerAdmin)


admin.site.register(Collection)

