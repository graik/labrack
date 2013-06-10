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
from labrack.models.generalmodels import Rack ,Container ,Location
from labrack.models.models import Attachement ,DnaComponent ,Cell ,CellSample, Unit

from labrack.models.models import Source ,DnaComponentType, PlasmidSample, CellType
from labrack.forms import DnaComponentForm, PlasmidSampleForm, CellForm,CellSampleForm

 


from labrack.models import utilLabrack
from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from django.contrib.sites.models import Site
 


class AutoAssignAdmin():

    # Save the owner of the object
    def save_model(self, request, obj, form, change):

        self.obj = obj  # store the obj for save_related method
        self.obj.createdBy = request.user
        #self.obj.owners.add(request.user)
        self.obj.save()


class userFilter(SimpleListFilter):
    title = 'User'
    parameter_name = 'user'

    def lookups(self, request, model_admin):
        listUser = User.objects.all()    
        r = [(c.id, c) for c in listUser]
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset
        else:
            return queryset.filter(createdBy__id__exact=self.value())
        
class ContainerAdmin(AutoAssignAdmin, admin.ModelAdmin):

    actions = ['make_csv']

    readonly_fields = ('registrationDate',)


    exportFields = OrderedDict( [('Container ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Container Type','containerType'),
                                 ('Rack','rack'),
                                 ('Created by','createdBy'),
                                 ('Description','description'),
                                 ('Registered at','creationDate'),
                                 ('last modified','modificationDate'),
                                 ])

    fieldsets = [
        (None, {
            'fields' : (('displayId', 'name'),
                        ('containerType', 'rack'),
                        'description',('createdBy','registrationDate')                        
                        ),
        }
         ),
    ]


    list_display = ('displayId', 'name', 'containerType', 'rack_url','location_url',
                    'createdBy')

    list_filter = ('containerType', 'rack__currentLocation__displayId', 
                   'rack__currentLocation__room', 
                   'rack__currentLocation__temperature', 'createdBy')

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
        if obj.rack.currentLocation:
            url = obj.rack.currentLocation.get_relative_url()
        else:
            url = ''

        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.rack.currentLocation))

    location_url.allow_tags = True

    location_url.short_description = 'Location'     


    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Container')


    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'createdBy':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(ContainerAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)

    #def get_readonly_fields(self, request, obj=None):
        #if obj is not None:
            #return self.readonly_fields + ('created_by',)
        #return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['createdBy'] = request.user
        request.GET = data
        return super(ContainerAdmin, self).add_view(request, form_url="", extra_context=extra_context)


class LocationAdmin(AutoAssignAdmin,admin.ModelAdmin):

    actions = ['make_csv']

    readonly_fields = ('registrationDate',)

    exportFields = OrderedDict( [('Location ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Temperature','temperature'),
                                 ('Room','room'),
                                 ('Registered at','creationDate'),
                                 ('Last modified','modificationDate'),
                                 ])

    fields = (('displayId', 'name'),'temperature','room','description',('createdBy','registrationDate'))

    list_display = ('displayId', 'name', 'temperature', 'room')

    list_filter = ('room', 'temperature')

    save_as = True

    search_fields = ('displayId', 'name',)

    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Location')

    make_csv.short_description = 'Export as CSV'

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'createdBy':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(LocationAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)


    #def get_readonly_fields(self, request, obj=None):
        #if obj is not None:
            #return self.readonly_fields + ('createdBy',)
        #return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['createdBy'] = request.user
        request.GET = data
        return super(LocationAdmin, self).add_view(request, form_url="", extra_context=extra_context)


class RackAdmin(AutoAssignAdmin,admin.ModelAdmin):

    actions = ['make_csv']

    readonly_fields = ('registrationDate',)


    exportFields = OrderedDict( [('Rack ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Current location','location_url'),
                                 ])

    list_display = ('displayId', 'name', 'location_url', 
                    )    
    fields = (('displayId', 'name'),'currentLocation','description',('createdBy','registrationDate'))

    def container_url(self, obj):
        url = obj.container.get_relative_url()
        return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.container.__unicode__()))

    def location_url(self, obj):
        if obj.currentLocation:
            url = obj.currentLocation.get_relative_url()
            locatAttr = obj.currentLocation.__unicode__()
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
        if db_field.name == 'createdBy':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(RackAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)

    #def get_readonly_fields(self, request, obj=None):
        #if obj is not None:
            #return self.readonly_fields + ('createdBy',)
        #return self.readonly_fields

    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['createdBy'] = request.user
        request.GET = data
        return super(RackAdmin, self).add_view(request, form_url="", extra_context=extra_context)


class categorieFilterForDnaComponentType(SimpleListFilter):
    title = 'Categorie'
    parameter_name = 'categorie'

    def lookups(self, request, model_admin):
        dnaType = DnaComponentType.objects.filter(subTypeOf=None)        
        r = [(c.name,c.name) for c in dnaType]
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset.exclude(subTypeOf=None)  
        else:
            return queryset.filter(subTypeOf__name=val)
    
        
class DnaComponentTypeAdmin(admin.ModelAdmin):

    fieldsets = (
        (None, {
            'fields': (('uri')
                ,('name',),
                       ('subTypeOf',)
                       ),
        }
         ),
    )

    list_display = ('name', 'categorie','categorieSubType',)
    list_filter = (categorieFilterForDnaComponentType,)
    ordering = ('subTypeOf',)
    

    search_fields = ('name', 'subTypeOf')

    def categorie(self, obj):
        if (obj.subTypeOf==None):
            categorieName = ''
        else:
            categorieName = obj.subTypeOf
        return mark_safe('%s' % (categorieName))
    categorie.allow_tags = True    
    categorie.short_description = 'Category'

    def categorieSubType(self, obj):
        return mark_safe('%s / %s' % (obj.name,obj.subTypeOf))
    categorieSubType.allow_tags = True    
    categorieSubType.short_description = 'Category / Type'

class ComponentAdmin(admin.ModelAdmin):


    actions = ['make_csv']

    exportFields = OrderedDict( [('ID', 'displayId'),
                                 ('Name', 'name'),
                                 ('Description','description'),
                                 ('Status','status'),
                                 ('Created by','createdBy'),
                                 ('registered at','creationDate'),
                                 ('modified at','modificationDate'),
                                 ])

    fieldsets = (
        (None, {
            'fields': (('displayId', 'name','status',))
        }
         ),
        ('Details', {
            'fields' : (('comment',),
                        ('attachements',)),            

        }
         ),
    )

    list_display = ('displayId', 'name', 'createdBy', 'showComment','status')
    list_filter = ('status',)

    ordering = ('displayId',)

    search_fields = ('displayId', 'name', 'comment')

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
        dnaOfMarker = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Marker')        
        r = [(c.id, c.displayId) for c in dnaOfMarker]
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset
        else:
            return queryset.filter(marker__id__exact=self.value())
 

class categorieFilterForPlasmid(SimpleListFilter):
    title = 'Categorie'
    parameter_name = 'categorie'

    def lookups(self, request, model_admin):
        dnaType = DnaComponentType.objects.filter(subTypeOf=None)        
        r = [(c.name,c.name) for c in dnaType]
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset
        else:
            dnaOfType = DnaComponent.objects.filter(componentSubType__subTypeOf__name=val)        
            return queryset.filter(dnaComponent__in=dnaOfType)
        
class categorieFilterForDna(SimpleListFilter):
    title = 'Categorie'
    parameter_name = 'categorie'

    def lookups(self, request, model_admin):
        dnaType = DnaComponentType.objects.filter(subTypeOf=None)        
        r = [(c.name,c.name) for c in dnaType]
        r.append((-1,'All including annotations'))
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == '-1'):
                        #return queryset.exclude(componentType__name__exact='Annotation')
            return queryset
        else:        
            if (val == None):
                #return queryset
                return queryset.filter(componentSubType__subTypeOf__name='Plasmid')
            else:
                return queryset.filter(componentSubType__subTypeOf__name=val)

class PlasmidSampleAdmin(admin.ModelAdmin):
    form = PlasmidSampleForm     

    readonly_fields = ('registrationDate',)

    actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']
    date_hierarchy = 'preparationDate'
    exportFields = OrderedDict([('Sample ID', 'displayId'),
                                 ('Location', 'container.location'),
                                 ('Container', 'container'),
                                 ('Number of aliquots', 'aliquotNr'),
                                 ('Status', 'status'),
                                 ('Comment', 'comment'),
                                 ('History', 'strProvenance()'),
                                 ('Sample link','sampleLinkStr()'),
                                 ('Preparation date', 'preparationDate'),
                                 ('Registered at','creationDate'),
                                 ('Last modified','modificationDate'),
                                 ])
    fieldsets = [
        (None, {
            'fields' : ((('container', 'displayId', 'status'),
                         ('preparationDate','createdBy','registrationDate',),
                         ('concentration','concentrationUnit','amount','amountUnit',),
                         ('solvent','aliquotNr',),
                         ('comment'),
                         )
                        )
        }
         ), 
        ('DNA Content',{'fields':[('dnaComponent','inCell'),]
                        }),
          
    ]
    list_display   = ('diplayId_readOnly', 'location_url', 'preparationDate',
                      'createdBy_User', 'content_url', 'concentraunit', 'amtunit',
                       'shorten_description','status','plasmidSample_editOnly')
    
    ordering       = ('container', 'displayId')

    save_as        = True

    save_on_top = True

    search_fields  = ('diplayId', 'name', 'container__displayId', 
                      'container__rack__currentLocation__displayId', 
                      'container__rack__currentLocation__location__temperature', 
                      'container__rack__currentLocation__location__room')

    list_filter = ComponentAdmin.list_filter.__add__((categorieFilterForPlasmid,))
    ##raw_id_fields = ComponentAdmin.raw_id_fields.__add__(('dnaComponent',))
    
    def diplayId_readOnly(self, obj):
        return mark_safe('<a href="%s/%s/%s/">%s</a>' % (S.admin_root,'../../reviewplasmid', obj.id,obj.displayId))
    diplayId_readOnly.allow_tags = True    
    diplayId_readOnly.short_description = 'ID'
 
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
    concentraunit.short_description = 'Concentration' 
    
    def amtunit(self,obj):
        if obj.amount == None:
            amount = ''
        else:
            amount = str(obj.amount)
            
        if obj.amountUnit == None:
            amountUnit = ''
        else:
            amountUnit = str(obj.amountUnit)        
        
        return amount + ' '+ amountUnit
    amtunit.short_description = 'Amount' 
    
    def location_url(self, obj):
            con = obj.container
            rac = con.rack
            locat = rac.currentLocation
            if locat:
                url = locat.get_relative_url()
                locatAttr = obj.container.rack.currentLocation.__unicode__()
            else:
                url = ''
                locatAttr = '' 
            return mark_safe('<a href="%s/%s">%s, %s</a>' % (S.admin_root, url, obj.container.rack.__unicode__(),locatAttr))
    location_url.allow_tags = True
    location_url.short_description = 'Location'  
    
    def shorten_description(self,obj):
        if not obj.comment:
            return u''
        if len(obj.comment) < 40:
            return unicode(obj.comment)
        return unicode(obj.comment[:38] + '..')
    shorten_description.short_description = 'Comment' 

    def content_url(self, obj):
        ##url = obj.dnaComponent.get_relative_url()
        ##return mark_safe('<a href="%s/%s">%s</a>' % (S.admin_root, url, obj.dnaComponent.__unicode__()))
        url = obj.dnaComponent.get_absolute_url()
        return mark_safe('<a href="../%s">%s</a>' % (url, obj.dnaComponent.__unicode__()))
    content_url.allow_tags = True
    content_url.short_description = 'Content'
    
    def createdBy_User(self, obj):
        url = obj.dnaComponent.get_relative_url()
        return mark_safe(obj.createdBy.__unicode__())
    createdBy_User.allow_tags = True
    createdBy_User.short_description = 'By'        
    
    def add_view(self, request, form_url="", extra_context=None):
        data = request.GET.copy()
        data['createdBy'] = request.user

        request.GET = data
        return super(PlasmidSampleAdmin, self).add_view(request, form_url="", extra_context=extra_context)
    
    def plasmidSample_editOnly(self, obj):
        return mark_safe('<a href="%s/%s/%s"><img src="http://icons.iconarchive.com/icons/custom-icon-design/office/16/edit-icon.png"/></a>' % (S.admin_root,'plasmidsample', obj.id))
    plasmidSample_editOnly.allow_tags = True    
    plasmidSample_editOnly.short_description = 'Edit'     


class DnaComponentAdmin(ComponentAdmin):
    form = DnaComponentForm


    readonly_fields = ('registrationDate',)


    fieldsets = (
        (None, {
            'fields': (('displayId', 'name','status'),
                       ('componentType','componentSubType','circular',),
                       ('vectorBackbone','marker','insert' ),
                       ('createdBy','registrationDate','source')
                       )
        }
         ),
        ('Details', {
            'fields' : (('comment',),
                        ('sequence'),
                        ('attachements',)
                        )
        }
         ),            
    )


    list_display = ('diplayId_readOnly', 'name', 'created_by_user','insert_url','vectorBackbone_url', 'marker_url', 'shorten_description','status','dnaComponent_editOnly')

    search_fields = ComponentAdmin.search_fields.__add__(('sequence',))
    list_filter = ComponentAdmin.list_filter.__add__(( vectorBackboneFilter, markerFilter, categorieFilterForDna,userFilter))



    def shorten_description(self,obj):
        if not obj.comment:
            return u''
        if len(obj.comment) < 40:
            return unicode(obj.comment)
        return unicode(obj.comment[:38] + '..')
    shorten_description.short_description = 'Comment'        

    def diplayId_readOnly(self, obj):
        return mark_safe('<a href="%s/%s/%s/">%s</a>' % (S.admin_root,'../../reviewdna', obj.displayId,obj.displayId))

        #return mark_safe('<a href="%s/%s/%s">%s</a>' % (S.search_root,'dnacomponent/objects', obj.id,obj.displayId))
    diplayId_readOnly.allow_tags = True    
    diplayId_readOnly.short_description = 'ID'

    def created_by_user(self,obj):
        return obj.createdBy
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
        ln = ""
        for m in obj.marker.all():
            url = m.get_absolute_url()
            if (ln==""):
                ln += mark_safe('<a href="../../%s">%s</a>' % (url, m.name))
            else:
                ln += " ,"+mark_safe('<a href="../../%s">%s</a>' % (url, m.name))
        for m in obj.vectorBackbone.marker.all():
            url = m.get_absolute_url()
            if (ln==""):
                ln += mark_safe('<a href="../../%s">%s</a>' % (url, m.name))
            else:
                ln += " ,"+mark_safe('<a href="../../%s">%s</a>' % (url, m.name))

        return ln        
        
    marker_url.allow_tags = True
    marker_url.short_description = 'Marker'   

    def insert_url(self, obj):
        if (obj.insert==None):
            return mark_safe('')
        else:
            url = obj.insert.get_absolute_url()
            return mark_safe('<a href="../../%s">%s</a>' % (url, obj.insert.name))
    insert_url.allow_tags = True
    insert_url.short_description = 'Insert'

    def dnaComponent_editOnly(self, obj):

        return mark_safe('<a href="%s/%s/%s"><img src="http://icons.iconarchive.com/icons/custom-icon-design/office/16/edit-icon.png"/></a>' % (S.admin_root,'dnacomponent', obj.id))

    dnaComponent_editOnly.allow_tags = True    
    dnaComponent_editOnly.short_description = 'Edit'    

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'createdBy':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        if db_field.name == 'vectorBackbone':
            dnaOfBaseVector = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Vector Backbone')        

            kwargs['queryset'] = dnaOfBaseVector

        if db_field.name == 'insert':
            #dnaOfInsert = DnaComponent.objects.filter(componentSubType__name='Insert')
            dnaOfInsert = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Insert')
            kwargs['queryset'] = dnaOfInsert

        return super(DnaComponentAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)

    def formfield_for_manytomany(self, db_field, request, **kwargs):
        
            
        if db_field.name == "componentType":
            kwargs['queryset'] = DnaComponentType.objects.filter(subTypeOf=None)

        if db_field.name == 'marker':
            kwargs['queryset'] = DnaComponent.objects.filter(componentSubType__subTypeOf__name='Marker')        


        return super(DnaComponentAdmin, self).formfield_for_manytomany(
            db_field, request, **kwargs)    


class categorieFilterForCellType(SimpleListFilter):
    title = 'Categorie'
    parameter_name = 'categorie'

    def lookups(self, request, model_admin):
        cellType = CellType.objects.filter(subTypeOf=None)        
        r = [(c.name,c.name) for c in cellType]
        return r

    def queryset(self, request, queryset):
        val = self.value()
        if (val == None):
            return queryset.exclude(subTypeOf=None)  
        else:
            return queryset.filter(subTypeOf__name=val)
        
class CellTypeAdmin(ComponentAdmin):
    fieldsets = (
            (None, {
                'fields': (('uri'),
                                  ('name'),('subTypeOf'),
                            )
            }
             ),           
    )

    list_display = ('name', 'categorie','categorieSubType',)
    list_filter = (categorieFilterForCellType,)
    
    ordering = ('subTypeOf',)

    search_fields = ('name', 'subTypeOf')

    def categorie(self, obj):
        if (obj.subTypeOf==None):
            categorieName = ''
        else:
            categorieName = obj.subTypeOf
        return mark_safe('%s' % (categorieName))
    categorie.allow_tags = True    
    categorie.short_description = 'Category'

    def categorieSubType(self, obj):
        return mark_safe('%s / %s' % (obj.name,obj.subTypeOf))
    categorieSubType.allow_tags = True    
    categorieSubType.short_description = 'Category / Type'
    
        
class CellAdmin(ComponentAdmin):
    form = CellForm
    
    actions = ['make_csv']
    
    readonly_fields = ('registrationDate',)
   
    fieldsets = (
            (None, {
                'fields': (('displayId', 'name','status'),
                           ('componentType','componentSubType',),
                           ('createdBy','registrationDate')
                            )
            }
             ),
            ('Details', {
                        'fields' : (('comment',),
                                    ('attachements',)
                                    )
                    }
                    ),             
    )
    
    list_display = ( 'name', 'createdBy', 'shorten_description','status')
    
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'createdBy':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(ChassisAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    def shorten_description(self,obj):
            if not obj.comment:
                return u''
            if len(obj.comment) < 40:
                return unicode(obj.comment)
            return unicode(obj.comment[:38] + '..')
    shorten_description.short_description = 'Comment'     

    #def get_readonly_fields(self, request, obj=None):
        #if obj is not None:
            #return self.readonly_fields + ('createdBy',)
        #return self.readonly_fields

    def make_csv(self, request, queryset):
        return importexport.generate_csv(self, request, queryset, 
                                         self.exportFields, 'Cell')
    make_csv.short_description = 'Export as CSV'    
    
    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == 'createdBy':
            kwargs['queryset'] = User.objects.filter(username=request.user.username)
        return super(CellAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)
    
    
class CellSampleAdmin(admin.ModelAdmin):
    form = CellSampleForm     

    readonly_fields = ('registrationDate',)

    actions = ['make_csv', 'make_ok', 'make_empty', 'make_bad']
    date_hierarchy = 'preparationDate'
    exportFields = OrderedDict([('Sample ID', 'displayId'),
                                 ('Location', 'container.location'),
                                 ('Container', 'container'), 
                                 ('Status', 'status'),
                                 ('Comment', 'comment'),  
                                 ('Preparation date', 'preparationDate'),
                                 ('Registered at','creationDate'), 
                                 ])
    fieldsets = [
        (None, {
            'fields' : ((('container', 'displayId', 'status'),
                         ('preparationDate','createdBy','registrationDate',),
                         ('concentration','concentrationUnit','amount','amountUnit',),
                         ('solvent','aliquotNr',),
                         ('comment'),
                         )
                        )
        }
         ), 
        ('DNA Content',{'fields':[('dnaComponent','inCell'),]
                        }),
          
    ]
    list_display   = ('diplayId_readOnly', 'location_url', 'preparationDate',
                      'createdBy_User', 'content_url', 'concentraunit', 'amtunit',
                       'shorten_description','status','plasmidSample_editOnly')
    
    ordering       = ('container', 'displayId')

    save_as        = True

    save_on_top = True

    search_fields  = ('diplayId', 'name', 'container__displayId', 
                      'container__rack__currentLocation__displayId', 
                      'container__rack__currentLocation__location__temperature', 
                      'container__rack__currentLocation__location__room')

    list_filter = ComponentAdmin.list_filter.__add__((categorieFilterForPlasmid,))

## Register regular admin panels
admin.site.register(Container, ContainerAdmin)
admin.site.register(Location, LocationAdmin)
admin.site.register(Rack, RackAdmin)
admin.site.register(DnaComponent, DnaComponentAdmin)
admin.site.register(DnaComponentType,DnaComponentTypeAdmin)
admin.site.register(Attachement)
admin.site.register(Source)
admin.site.register(PlasmidSample,PlasmidSampleAdmin)
admin.site.register(Cell,CellAdmin)
admin.site.register(CellSample)
admin.site.register(Unit)

admin.site.register(CellType,CellTypeAdmin)






