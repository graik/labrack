from usermixin import UserMixinModel
from django.db import models
from datetime import datetime
from django.utils.safestring import mark_safe

from labrack.models.generalmodels import Rack
from labrack.models.generalmodels import Container
from labrack.models.generalmodels import Location
from labrack.models.generalmodels import Unit

from django.contrib.auth.models import User


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

    comment = models.TextField( 'Detailed description', blank=True)
    
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
        if not self.comment:
            return u''
        if len(self.comment) < 40:
            return unicode(self.comment)
        return unicode(self.comment[:38] + '..')
    showComment.short_description = 'Comment'

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

    name = models.CharField('Name', unique=True, max_length=200, 
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
        link = self.uploadedFile if self.uploadedFile else ''
        name = link
        if str(name)<>"":
            names=str(name).split("/")
            name = names[len(names)-1]
        return u'%s' % (link)
    
    
    
    def get_absolute_html_url(self):
        link = self.uploadedFile if self.uploadedFile else ''
        name = link
        if str(name)<>"":
            names=str(name).split("/")
            name = names[len(names)-1]
        return u'<a href="../../media/%s">%s</a>' % (link,name)         
    
    class Meta:
        app_label = 'labrack' 
        verbose_name = 'Extra Files'    
    

class DnaComponent(Component):
    """
    Description of a stretch of DNA or RNA.
    """
   
    sequence = models.TextField( help_text='nucleotide sequence', blank=True, null=True )      
    
    componentSubType = models.ManyToManyField(DnaComponentType, 
                                           verbose_name='Type',
                                           related_name='Type',  blank=True, null=True)    
    
    circular = models.BooleanField( 'Circular', default=True, 
                                    help_text='is the DNA Circular or not.')
    attachements = models.ManyToManyField(ExtraFile,blank=True,related_name='Attachements', null=True)
    
    source = models.ForeignKey( Source, blank=True, null=True, 
                                      help_text='external source' )
    
    vectorBackbone = models.ForeignKey( 'self', blank=True, null=True ,verbose_name='Vector Backbone')
    
    marker = models.ManyToManyField( 'self', blank=True, null=True, 
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
    
    comment = models.TextField('Comments', blank=True)



    preparation_date = models.DateField(default=datetime.now(), verbose_name="created")
    
    registration_date = models.DateField(default=datetime.now(), verbose_name="registred")
    
    historyDescription = models.TextField('Comment', null=True,blank=True,max_length = 255) 
    
    attachements = models.ManyToManyField(ExtraFile,blank=True, null=True)
    
    source = models.ForeignKey( Source, blank=True, null=True, 
                                      help_text='external source' )
    
    solvent = models.CharField('Buffer/Medium', max_length=100, blank=True)

    concentration = models.FloatField('Concentration', null=True, blank=True)

    concentrationUnit = models.ForeignKey(Unit, 
                                          verbose_name='Concentration Unit',
                                          related_name='concUnit',
                                          null=True, blank=True,
                                          limit_choices_to = {'unitType': 'concentration'})    
    
    amount = models.FloatField('Amount', null=True, blank=True)
    amountUnit = models.ForeignKey(Unit, 
                                          verbose_name='Amount Unit',
                                          related_name='amouUnit', 
                                          null=True, blank=True, 
                                          limit_choices_to = {'unitType': 'concentration'})
    
    
    
    

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
        
        
### To be removed what under

class Background(models.Model):
    
    displayId = models.CharField('displayId', max_length=100, unique=True, 
        help_text='Unique identification')
    description = models.TextField('Comments', blank=True)
    class Meta:
        app_label = 'labrack'
        
class Yeast(models.Model):
    
    
    mtid = models.CharField('MTID', max_length=20, unique=True, 
        help_text='Unique identification')
    
    
    MATING_TYPE_CHOICES = (('a','a'),('alpha','alpha'),('a/alpha','a/alpha'))
    matingType = models.CharField('matingType',max_length=200, null=True, blank=True, choices=MATING_TYPE_CHOICES)
    
    
    background = models.ForeignKey(Background,
                            verbose_name='Background',
                            null=True, blank=True)
    
    GENO_CHOICES = (('ade2','ade2'),('his3','his3'),('leu2','leu2'),('lys2','lys2'),('trp1','trp1'),('ura3','ura3'),('gal','gal'),('kan','kan'),('nat','nat'))
    genoList = models.CharField('genoList',max_length=200, null=True, blank=True, choices=GENO_CHOICES)
    
   
    
    relevantGeno = models.TextField('RelevantGeno', null=True,blank=True) 
    fillingDate = models.DateField(default=datetime.now(), null=True,blank=True)
    user = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='%(class)s_created_by')
    
    source_old_list = models.ForeignKey(Source,
                            verbose_name='source',
                            null=True, blank=True)
    
    sourceLabMember = models.ForeignKey(User, null=True, blank=True)
    
    attachedSheet = models.BooleanField('Attached sheet',default=False)
    
    otherMarkers = models.TextField('otherMarkers', null=True,blank=True)
    
    strainConstruction = models.TextField('strainConstruction', null=True,blank=True)
    
    comments = models.TextField('Comments', null=True,blank=True) 
     
    fileBrowse1 = models.FileField(upload_to='documents/%Y/%m/%d',blank=True,null=True)    
    fileBrowse2 = models.FileField(upload_to='documents/%Y/%m/%d',blank=True,null=True)
    
    
    def __unicode__(self):
        return self.mtid
    class Meta:
        app_label = 'labrack'
