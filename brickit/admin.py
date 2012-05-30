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

from brickit.models import *
from django.contrib import admin
from django.http import HttpResponse

class ComponentAdmin(admin.ModelAdmin):
    """
    Abstract base class with common description for DNA / protein / cell
    forms
    """
    save_as        = True # activate 'Save As' button

    fieldsets = (
        (None,
         {'fields': (('displayId', 'name', 'uri'),
                     ('shortDescription','abstract'),
                     ('users',))}),
        ('Details',
         {'fields' : ('componentType','variantOf','description',)}),
        ('Annotations',
         {'fields': ('annotations'),
          'classes':('collapse',)} ),
        )

    list_display   = ('displayId','name','componentType', 'shortDescription',
##                      'isavailable_cells',
##                      'isavailable_dna'
                      )
    list_filter    = ('users', 'status', 
                      'componentType',
                      )
    ordering       = ('displayId',)
    search_fields  = ('displayId','name','shortDescription','description',
                      'users__username' )

class DnaComponentAdmin(admin.ModelAdmin):

    fieldsets = (
        (None,
         {'fields': (('displayId', 'name', 'uri'),
                     ('shortDescription','abstract'),
                     ('users','status'))}),
        ('Details',
         {'fields' : (('componentType','variantOf'),
                     ('optimizedFor', 'translatesTo'),
                     ('description',))}),
        ('Sequence Details',
         {'fields': ('sequence','annotations'),
          'classes':('collapse',)} ),
        )

    list_display   = ('displayId','name','shortDescription',
                      'optimizedFor','translatesTo', 'status'
##                      'isavailable_cells',
##                      'isavailable_dna'
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


class SampleAdmin(admin.ModelAdmin):
    fieldsets = (
        (None, {
        'fields' : (('displayId',), 
                    ('container','vesselType','users'),
                    ('sampleType',),
                    ('concentration', 'concentrationUnit'),
                    ('comments',)
                    )
        }),
    )

    list_display   = ('displayId', 'created', 'container',
                      'concentration' )
    list_filter    = ('vesselType','sampleType','container','created',
                      'users')
    ordering       = ('displayId',)
    search_fields  = ('displayId','comments','users__username',
                      'container__displayId')
    
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

