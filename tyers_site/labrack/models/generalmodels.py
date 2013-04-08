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

from permissionModel import PermissionModel



## importing custom models

import os


from datetime import datetime

import urllib
import tyers_site.settings as S


class Location(PermissionModel):
    """
    A location (fridge, freezer, shelf) where containers are stored
    """
    displayId = models.CharField('Location ID', max_length=20, unique=True, 
                                 help_text='Unique identifier')

    name = models.CharField('Name', max_length=200, 
                            blank=True, help_text='Informative name of fridge, freezer or shelf')

    temperature = models.FloatField('Temperature', help_text= unichr(176) +'C')

    room = models.CharField('Room', max_length=20)

    #creation_date = models.DateTimeField('Created at', auto_now_add=True)

    #modification_date = models.DateTimeField('Modified at', auto_now=True)

    description = models.CharField( 'short Description', max_length=200,
                                            help_text='Very short description for listings')    

    registration_date = models.DateField(default=datetime.now(), verbose_name="registred")
   

    def __unicode__(self):
        return self.displayId + " (" + self.temperature.__str__() \
               + unichr(176) + "C)  [Room: " + self.room + "]"

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'location/%i/' % self.id


    def related_racks(self):
        r = Rack.objects.filter( current_location = self.id )
#        s = [content.sample for content in r]        

        return r 

    class Meta:
        app_label = 'labrack'
        
          
        
        





class Rack(PermissionModel):
    """
    A Rack (box) where containers are stored
    """
    displayId = models.CharField('Rack ID', max_length=20, unique=True, 
                                 help_text='Unique identifier')

    name = models.CharField('Name', max_length=200, 
                            blank=True, help_text='Informative name of rack ex : R1-F1 for a rack named R1 that will not move from Fridge 1')

    current_location = models.ForeignKey( 'Location', blank=True, null=True )
    
    description = models.CharField( 'short Description', max_length=200,
                                            help_text='Very short description for listings')    
    
    registration_date = models.DateField(default=datetime.now(), verbose_name="registred")
       
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

    displayId = models.CharField('Box ID', max_length=20, unique=True, 
                                 help_text='Unique identifier. E.g. D001 or C012 (D...DNA, C...Cells)')

    name = models.CharField('Name', max_length=200, blank=True, 
                            help_text='Informative name for tables and listings')

    containerType = models.CharField('Type of container', max_length=30, 
                                     choices=STORAGE_CONTAINER_TYPES )

    rack = models.ForeignKey(Rack)

    registration_date = models.DateField(default=datetime.now(), verbose_name="registred")
       
    #: optional long description
    description = models.TextField( 'Detailed description', blank=True)

    #creation_date = models.DateTimeField('Created at', auto_now_add=True)

    #modification_date = models.DateTimeField('Modified at', auto_now=True)

    # rack = models.ForeignKey(Rack, null=True, blank=True)

    def __unicode__( self ):
        return "%s  ,%s, %s" % (self.displayId,self.rack,self.rack.current_location)

    def related_dna_samples( self ):
        """
        """

        import labrack.models as M
        
        r = M.DnaSample.objects.filter(container=self.id)
#        s = [content.sample for content in r]        

        return r
    
    def related_chassis_samples( self ):
        """
        """

        import labrack.models as M
        
        r = M.ChassisSample.objects.filter(container=self.id)
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

    displayId = models.CharField('ID', max_length=20, unique=True, 
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

   
    #precedes = models.ManyToManyField( 'self', symmetrical=False,
    #                                   null=True, blank=True )

   
    # override Django object methods
    def __unicode__(self):
        return u'(%s - %s) (%s,%s)' % (self.subComponent, self.componentAnnotated,self.bioStart,self.bioEnd)

    

    class Meta:
        app_label = 'labrack'    
        abstract = True

class DnaSequenceAnnotation(SequenceAnnotation):
    subComponent = models.ForeignKey('DnaComponent', related_name ='subComponentOf')

    componentAnnotated = models.ForeignKey('DnaComponent', related_name ='annotatedForComponent',null=True, blank=True)
    
    strand = models.CharField('Strand',max_length=1, blank=True, null=True,
                              choices=(('+','+'),('-','-')) )
    
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'dnasequenceannotation/%i/' % self.id    
    class Meta:
        app_label = 'labrack'    
    
    @staticmethod
    def isRelated(displayId1,displayId2):
        from labrack.models.component import DnaComponent
        try:
            id1 = DnaComponent.objects.get(displayId=displayId1)
            id2 = DnaComponent.objects.get(displayId=displayId2)
        
            r1 = DnaSequenceAnnotation.objects.filter( subComponent=id1,componentAnnotated=id2).count()
            r2 = DnaSequenceAnnotation.objects.filter( subComponent=id2,componentAnnotated=id1).count()
        
            if (r1>0 or r2>0):
                return True
        except error:
            return False
        return False
    
    @staticmethod
    def deleteRelated(displayId):
        from labrack.models.component import DnaComponent
        
        id = DnaComponent.objects.get(displayId=displayId)                
        DnaSequenceAnnotation.objects.filter( subComponent=id).delete()
        
        
        


class ProteinSequenceAnnotation(SequenceAnnotation):
    subComponent = models.ForeignKey('ProteinComponent', related_name ='subComponentOf')

    componentAnnotated = models.ForeignKey('ProteinComponent', related_name ='annotatedForComponent',null=True, blank=True)
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'proteinsequenceannotation/%i/' % self.id    
    class Meta:
        app_label = 'labrack'    
