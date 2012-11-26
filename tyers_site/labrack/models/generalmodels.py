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


from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
from django.core.files.storage import FileSystemStorage
from django.contrib import messages
from django.db import models
from django.db.models import Q
from django.utils.safestring import mark_safe
from django.utils.safestring import SafeUnicode
from Bio import SeqIO
from Bio.Seq import Seq
from tyers_site import settings   

from tyers_site.labrack.models.component import *
from permissionModel import PermissionModel







## importing custom models

import os


from datetime import datetime

import urllib
import tyers_site.settings as S


class Location(models.Model):
    """
    A location (fridge, freezer, shelf) where containers are stored
    """
    displayId = models.CharField(max_length=20, unique=True, 
                                 help_text='Unique identifier')

    name = models.CharField('Name', max_length=200, 
                            blank=True, help_text='Informative name of fridge, freezer or shelf')

    temperature = models.FloatField('Temperature', help_text= unichr(176) +'C')

    room = models.CharField('Room', max_length=20)

    creation_date = models.DateTimeField('Created at', auto_now_add=True)

    modification_date = models.DateTimeField('Modified at', auto_now=True)


    #rack = models.CharField('Rack', max_length=200, 
    #    blank=True, help_text='Informative name of rack')


    def __unicode__(self):
        return self.displayId + " (" + self.temperature.__str__() \
               + unichr(176) + "C)  [Room: " + self.room + "]"

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'location/%i/' % self.id




    class Meta:
        app_label = 'labrack'





class Rack(models.Model):
    """
    A Rack (box) where containers are stored
    """
    displayId = models.CharField(max_length=20, unique=True, 
                                 help_text='Unique identifier')

    name = models.CharField('Name', max_length=200, 
                            blank=True, help_text='Informative name of rack ex : R1-F1 for a rack named R1 that will not move from Fridge 1')

    current_location = models.ForeignKey( 'Location', blank=True, null=True )

    def __unicode__(self):
        return u'%s' % (self.displayId)

    def related_containers( self ):

        r = Container.objects.filter( rack = self.id )
#        s = [content.sample for content in r]        

        return r    

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'rack/%i/' % self.id          
    class Meta:
        app_label = 'labrack'   


class Container( PermissionModel ):
    """
    A container holding several physical samples of nucleic acids, proteins 
    or other stuff.
    """

    STORAGE_CONTAINER_TYPES= (
        ('96-well-plate', '96 well plate'),
        ('384-well-plate','384 well plate'),
        ('box', 'Freezer box'),
        ('other', 'other' ) )

    displayId = models.CharField(max_length=20, unique=True, 
                                 help_text='Unique identifier. E.g. D001 or C012 (D...DNA, C...Cells)')

    name = models.CharField('Name', max_length=200, blank=True, 
                            help_text='Informative name for tables and listings')

    containerType = models.CharField('Type of container', max_length=30, 
                                     choices=STORAGE_CONTAINER_TYPES )

    rack = models.ForeignKey(Rack)

    #: optional long description
    description = models.TextField( 'Detailed description', blank=True)

    creation_date = models.DateTimeField('Created at', auto_now_add=True)

    modification_date = models.DateTimeField('Modified at', auto_now=True)

    # rack = models.ForeignKey(Rack, null=True, blank=True)

    def __unicode__( self ):
        return "%s" % self.displayId

    def related_samples( self ):
        """
        """


        r = Sample.objects.filter( container = self.id )
#        s = [content.sample for content in r]        

        return r

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'container/%i/' % self.id    

    def next_well( self ):
        """
        Try to guess next free well.
        """
        r = ''

        if self.container_type == 'box':
            try:
                wells = [ int(s.well) for s in self.samples.all() ]
                wells.sort()
                r = str( wells[-1] + 1 )
            except:
                r = "%02i" % self.samples.count() + 1

        if self.container_type in ['96-well-plate', '384-well-plate']:
            try:
                wells = [ s.well for s in self.samples.all() ]
                wells.sort()
                row = wells[-1][0]
                col = wells[-1][1:] + 1
                r = "%s%02i" % (row, col)
            except:
                r = "%02i" % self.samples.count() + 1

        return r


    class Meta:
        ordering = ('displayId',)

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'container/%i/' % self.id





    class Meta:
        app_label = 'labrack'   





class SelectiveMarker( models.Model ):
    """
    Describes an Antibiotic or similar resistence marker
    """

    displayId = models.CharField('Component (ID)', max_length=20, unique=True, 
                                 help_text='Unique identification')

    name = models.CharField('Name', max_length=200, blank=True, 
                            help_text='Descriptive name (e.g. AMP)')

    description = models.TextField( 'Detailed description', blank=True)

    def __unicode__(self):
        return self.name

##    def get_absolute_url(self):
##        """Define standard URL for object.get_absolute_url access in templates """
##        return APP_URL+'/selectivemarker/%i/' % self.id

    class Meta:
        app_label = 'labrack'
        ordering = ('name',)





#class VectorDnaComponent(NucleicAcidComponent):
#    """
#    Description of vector backbone. So far identical to NucleicAcidComponent.
#    """
#
#    #: [optional] link to resistence or similar marker
#    marker = models.ManyToManyField(SelectiveMarker,
#                                     verbose_name='selective marker(s)',
#                                     null=True, blank=True)
#
#    class Meta(Component.Meta):
#        pass
#    
#    def get_relative_url(self):
#        """
#        Define standard relative URL for object access in templates
#        """
#        return 'vectordnacomponent/%i/' % self.id
#
#    def related_samples(self):
#        """
#        """
#        return self.vector_samples.all()

class Document(models.Model):
    docfile = models.FileField(upload_to='documents/BulkSamples/%Y/%m/%d')
    class Meta:
        app_label = 'labrack'

class SequenceAnnotation(models.Model):
    """
    Identify features and sub-components on a DNA or protein sequence.
    """
    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=False,
                           help_text='Specify external URI or leave blank for local URI')

    bioStart = models.PositiveIntegerField( 'From position', blank=True, null=True,
                                            validators=[],
                                            help_text='starting position counting from 1')

    bioEnd = models.PositiveIntegerField( 'To position', blank=True, null=True,
                                          validators=[],
                                          help_text='End position counting from 1')

    strand = models.CharField('Strand',max_length=1, blank=True, null=True,
                              choices=(('+','+'),('-','-')) )

    precedes = models.ManyToManyField( 'self', symmetrical=False,
                                       null=True, blank=True )

    subComponent = models.ForeignKey('Component', related_name ='subComponentOf')

    componentAnnotated = models.ForeignKey('Component', related_name ='annotatedForComponent',null=True, blank=True)

    # override Django object methods
    def __unicode__(self):
        return u'(%s - %s)(%s,%s)' % (self.subComponent, self.componentAnnotated,self.bioStart,self.bioEnd)

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'sequenceannotation/%i/' % self.id

    class Meta:
        app_label = 'labrack'    

