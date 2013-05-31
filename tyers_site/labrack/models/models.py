from usermixin import UserMixinModel
from django.db import models
from datetime import datetime
from django.utils.safestring import mark_safe

from labrack.models.generalmodels import Rack
from labrack.models.generalmodels import Container
from labrack.models.generalmodels import Location


class Component(UserMixinModel,models.Model):
    """
    Base class for cells, nucleic acids, proteins, and chemicals.
    Not shown to the user (currently) but the table exists and collects
    all fields that are common to all types of Components.
    
    See Meta.abstract
    """

    STATUS_CHOICES = ( ('available', 'available'),
                       ('planning', 'planning'),
                       ('under_construction', 'under construction'),
                       ('abandoned', 'abandoned'))

    displayId = models.CharField('ID',  max_length=20, unique=True, 
        help_text='Unique identification')

    name = models.CharField('Name', max_length=200, blank=True, 
                            help_text='Descriptive name (e.g. "TEV protease")')

    description = models.TextField( 'Detailed description', blank=True)
    
    registration_date = models.DateField(default=datetime.now(), verbose_name="registred")
   

    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=False, 
                           help_text='Specify external URI or leave blank for local URI')

    #: non-SBOL compliant status    
    status = models.CharField( max_length=30, choices=STATUS_CHOICES, 
                               default='planning')


    def formatedUrl(self):
        name = self.name if self.name else ''
        return SafeUnicode("<a href='/admin/labrack/" +self.get_relative_url() \
                           + "'>" + self.displayId + ' - ' +  name + "</a>")
    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s - %s' % (self.displayId, name)


    def showComment( self ):
        """
        @return: str; truncated comment
        """
        if not self.description:
            return u''
        if len(self.description) < 40:
            return unicode(self.description)
        return unicode(self.description[:38] + '..')
    showComment.short_description = 'comment'

    class Meta:
        app_label = 'labrack'
        abstract = False

class ComponentType( models.Model ):
    """
    Helper class for classifying parts.
    Following SBOL, each type should in theory correspond to a Sequence Ontology term.
    """

    uri = models.URLField( unique=False, blank=True, null=True,
                           help_text='Typically a sequence ontology URI, example: http://purl.obolibrary.org/obo/SO_0000167' )

    name = models.CharField('Name', max_length=200, 
                            help_text='Informative name')

    def __unicode__( self ):
        return unicode(self.name)

    class Meta:
        app_label = 'labrack' 
        abstract = True

class DnaComponentType( ComponentType ):

    #: required ? directional relationship to parent type or types
    subTypeOf = models.ForeignKey('self', blank=True, 
                                       null=True, related_name='subTypes')
        
    class Meta:
        app_label = 'labrack' 
        verbose_name = 'DNA Type'
        abstract = False

class Source(models.Model):
    name = models.CharField('Name', max_length=200, 
                            help_text='Informative name')
    description = models.TextField( 'Detailed description', blank=True)
    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s' % (self.name)
    
    class Meta:
            app_label = 'labrack' 
    
    
    
class ExtraFile(models.Model):
     
    FILE_TYPE_CHOICES = ( ('gbk', 'gbk'),
                           ('attachement1', 'attachement1'),
                           ('attachement2', 'attachement2'),
                           ('others', 'others'))
    uploadedFile = models.FileField(upload_to='documents/%Y/%m/%d',blank=True,null=True)

    fileType = models.CharField( max_length=30, choices=FILE_TYPE_CHOICES, 
                                   default='Type')    
    
    description = models.TextField( 'Detailed description', blank=True)
    
    def __unicode__(self):
        name = self.uploadedFile if self.uploadedFile else ''
        return u'%s' % (self.uploadedFile)
          
    
    class Meta:
        app_label = 'labrack' 
        verbose_name = 'Extra Files'    
    

class DnaComponent(Component):
    """
    Description of a stretch of DNA or RNA.
    """
   
    sequence = models.TextField( help_text='nucleotide sequence', blank=True, null=True )      
    
    
    #componentType = models.ManyToManyField(DnaComponentType, verbose_name='Categorie')        

    componentSubType = models.ManyToManyField(DnaComponentType, 
                                           verbose_name='Type',
                                           related_name='Type',  blank=True, null=True)    
    
    circular = models.BooleanField( 'Circular', default=True, 
                                    help_text='is the DNA Circular or not.')
    attachements = models.ManyToManyField(ExtraFile,blank=True,related_name='Attachements', null=True)
    
    source = models.ForeignKey( Source, blank=True, null=True, 
                                      help_text='external source' )
    
    baseVector = models.ForeignKey( 'self', blank=True, null=True )
    
    marker = models.ForeignKey( 'self', blank=True, null=True, 
                                      related_name='Marker')
    
    insert = models.ForeignKey( 'self', blank=True, null=True, 
                                      related_name='Insert')
    

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'dnacomponent/%i/' % self.id
    
    def get_absolute_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return '../../reviewdna/%s/' % self.displayId

          
    class Meta:
        app_label = 'labrack'
        verbose_name = 'DNA'


    
class Sample( UserMixinModel ,models.Model):
    

    displayId = models.CharField('ID', max_length=20,
                                 help_text='Label or well position. Must be unique within container.')

    name = models.CharField('Name', max_length=200, blank=True, 
                            help_text='Informative name for tables and listings')

    #: link to a single container
    container = models.ForeignKey( Container, related_name='samples' 
                                   )

    aliquotNr = models.PositiveIntegerField('Number of aliquots', 
                                            null=True, blank=True)

    STATUS_CHOICES = (('ok', 'ok'), 
                      ('preparation', 'in preparation'),
                      ('empty', 'empty'),
                      ('bad', 'bad'),
                      )
    status = models.CharField( max_length=30, choices=STATUS_CHOICES, 
                               default='ok')
    reference_status = models.BooleanField('Reference sample',default=False,
                                           help_text='mark sample as master or reference')


    description = models.TextField('Description / Comments', blank=True)



    preparation_date = models.DateField(default=datetime.now(), verbose_name="created")
    
    registration_date = models.DateField(default=datetime.now(), verbose_name="registred")
    
    historyDescription = models.TextField('Comment', null=True,blank=True,max_length = 255) 
    
    attachements = models.ManyToManyField(ExtraFile,blank=True, null=True)
    
    source = models.ForeignKey( Source, blank=True, null=True, 
                                      help_text='external source' )
    

class PlasmidSample(Sample):
    """
    Description of a stretch of DNA or RNA.
    """
    hostGenotype = models.TextField('host Genotype', null=True,blank=True,max_length = 255) 
    
    parentalVector = models.TextField('host Genotype', null=True,blank=True,max_length = 255) 
    
    synonyms = models.TextField('synonyms', null=True,blank=True,max_length = 255) 
    
    attachedSheet = models.BooleanField('Attached Sheet',default=False,
                                           help_text='is there a sheet document referencing the sample')
    
    dnaComponent = models.ForeignKey( DnaComponent, related_name='dna')
    
    class Meta:
        app_label = 'labrack'
        verbose_name = 'PlasmidSample'    