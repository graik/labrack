from django.db import models

from django.contrib.auth.models import User

import django.core.validators as V

# Create your models here.
class SampleContainer( models.Model ):
    """
    A container holding several physical samples of DNA or other stuff.
    """

    STORAGE_CONTAINER_TYPES= (
        ('96-well-plate', '96 well plate'),
        ('384-well-plate','384 well plate'),
        ('box', 'freezer box'),
        ('other', 'other' ) )
    
    #: actual label stuck to this container
    displayId = models.CharField('Label (ID)', max_length=20, unique=True,
                                 help_text='example: D001 or C012 (D...DNA, C...Cells)')

    #: [optional] annotation
    shortDescription = models.CharField(max_length=200, blank=True )

    #: what type of container is it
    containerType = models.CharField('type of container', max_length=30,
                                      choices=STORAGE_CONTAINER_TYPES )

    #: creators or users
    users = models.ManyToManyField(User, null=True, blank=True, db_index=True)

    #: optional long description
    description = models.TextField( 'Detailed description', blank=True,
        help_text='Use restructured text markup for links, lists and headlines')

    def __unicode__( self ):
        return "%s" % self.displayId

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


class Sample( models.Model ):
    """
    Sample describes a single tube or well holding DNA, cells or protein.
    """

    STORAGE_WELL_TYPES = ( ('tube','tube'), ('well','well'), ('other','other'))
    STORAGE_TYPES      = ( ('dna', 'DNA'), ('cells','cells'), ('protein','protein') )
    CONCENTRATION_UNITS= ( ('mg/l', 'mg/l'), ('mol/l','mol/l'), 
                           ('umol/l','umol/l') )

    #: actual label on the tube, enforced to be unique in combination with
    #: container
    displayId = models.CharField('ID (label)', max_length=20,
        help_text='label or well position. Must be unique within container.' )

    #: link to a single container
    container = models.ForeignKey( SampleContainer, related_name='samples' )

    wellType = models.CharField('type of well or tube', max_length=30,
                                 choices=STORAGE_WELL_TYPES,
                                 default='tube')

    sampleType = models.CharField('stored as', max_length=100,
                                    choices=STORAGE_TYPES, default='dna' )

    #: creators or users
    users = models.ManyToManyField(User, db_index=True)

    #: [optional]
    comments = models.TextField(blank=True)

    created = models.DateField('created at', auto_now_add=True)

##    #: link to the physical DNA contained in this sample
##    dna = models.ForeignKey( 'DNA',
##                             verbose_name='physical DNA',
##                             related_name='in_samples',
##                             blank=False)
##
##    #: [optional] strain this DNA is stored in, if any
##    cell = models.ForeignKey( 'Cell',
##                              verbose_name='in cell',
##                              related_name='sample',
##                              blank=True,
##                              null=True)

    #: concentration in ng / ul (=mg/l)
    concentration = models.DecimalField( max_digits=6, decimal_places=2,
                                         default=50 )
    #: unit of given concentration
    concentrationUnit = models.CharField( 'conc. unit', max_length=10,
                                           choices=CONCENTRATION_UNITS,
                                           default='mg/l')

    def __unicode__( self ):
        return self.container.displayId + ' / ' + self.displayId

##    def get_absolute_url(self):
##        """Define standard URL for object.get_absolute_url access in templates """
##        return APP_URL+'/sample/%i/' % self.id

##    def get_sequencing_evaluation( self ):
##        """
##        @return: str, best registered sequencing analysis 
##                ('confirmed', 'inconsistent', 'ambiguous' or 'not analyzed')
##        """
##        results = [ seq.evaluation for seq in self.sequencing.all() ]
##
##        for verdict in Sequencing.EVALUATION:
##            if verdict[0] in results:
##                return verdict[1]
##
##        return 'none'
##    get_sequencing_evaluation.short_description = 'Sequencing'

##    def related_samples( self ):
##        """
##        @return: QuerySet of samples with similar DNA content
##        """
##        if self.dna.biobrick:
##            r = Sample.objects.filter( dna__biobrick=self.dna.biobrick )
##        else:
##            r = Sample.objects.filter( dna__biobrick=self.dna.biobrick,
##                                       dna__vector=self.dna.vector )
##            
##        return r.exclude( pk=self.pk )
    
##    def get_sequence( self, recenter=0 ):
##        return self.dna.get_sequence( recenter=recenter )
##
##    def get_pretty_sequence( self, recenter=0 ):
##        return self.dna.get_pretty_sequence( recenter=recenter )

    class Meta:
        unique_together = ('displayId', 'container')
        ordering = ('container', 'displayId')


class ComponentType( models.Model ):
    """
    Helper class for classifying parts.
    Following SBOL, each type should correspond to a Sequence Ontology term.
    """
    
    #: required
    uri = models.URLField( unique=True,
        help_text='Typically a sequence ontology URI, example: '+\
        'http://purl.obolibrary.org/obo/SO_0000167' )
    
    #: required
    name = models.CharField( max_length=50,
                             help_text='(Sequence Ontology) name')

    #: required ? directional relationship to parent type or types
    subTypeOf = models.ManyToManyField('self', symmetrical=False,
                                       blank=True, null=True,
                                       related_name='subTypes')

    def __unicode__( self ):
        return self.name

class Component(models.Model):
    """
    Abstract base class for DNA and protein 'parts' as well as host cells.
    """
    
    #: required ID
    displayId = models.CharField('ID', max_length=20,
                    help_text='Identification; Unique within collection.', 
                    unique=True)
    
    #: optional short name
    name = models.CharField('Name', max_length=50, blank=True, null=True,
                            help_text='descriptive name (e.g. "TEV protease")')
    
    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=True,
                help_text='Specify external URI or leave blank for local URI')
    
    #: required short description -- will be added as first line in SBOL
    #: exports (where this field doesn't exist)
    shortDescription = models.CharField( 'short Description', max_length=200,
        help_text='very short description for listings')

    #: optional long description
    description = models.TextField( 'Detailed description', blank=True,
        help_text='Use restructured text markup for links, lists and headlines')
    
    #: optional type
    componentType = models.ManyToManyField( ComponentType, 
                                     blank=True, null=True,
                                     verbose_name='Type',
##                                     related_name='%(class)s',
                                     help_text='Classification of this part.')

    #: optional directional relationship to parent part or parts
    variantOf = models.ManyToManyField( 'self', symmetrical=False,
                                        blank=True, null=True,
                                        verbose_name='Variant of',
##                                        related_name='%(class)Variants',
        help_text='Specify part(s) this part is derived from, if any.' )
    
    #: optional classification as abstract part with missing information
    abstract = models.BooleanField( 'Abstract Part', default=False,
        help_text='Entry only serves as container to organize related parts.')
    
    #: creators or users authorized to change
    users = models.ManyToManyField(User, db_index=True,
                                   help_text='Primary users of this part.')
    
    annotations = models.ManyToManyField( 'SequenceAnnotation',
                                          verbose_name='Annotations',
                                          blank=True, null=True )
    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s %s' % (self.displayId, name)
    
    class Meta:
        abstract = True
    

class DnaComponent(Component):
    """
    Description of a stretch of DNA.
    """
    #: optional sequence
    sequence = models.TextField( help_text='nucleotide sequence', 
                                 blank=True, null=True )
    
    translatesTo = models.ForeignKey( 'ProteinComponent', 
                                      blank=True, null=True,
                                      related_name='encodedBy',
        help_text='Protein part this sequence translates to' )
    
    optimizedFor = models.ForeignKey( 'Chassis',
                                     blank=True, null=True )

    class Meta(Component.Meta):
        pass
    
    
class ProteinComponent(Component):
    """
    Description of a amino-acid encoded protein 'part'.
    """
    #: optional sequence
    sequence = models.TextField( help_text='amino acid sequence', 
                                 blank=True, null=True )

    class Meta(Component.Meta):
        pass


class Chassis(Component):
    """
    Description of a host system. Usually this will be a cell type or bacterial
    strain.
    """

    class Meta(Component.Meta):
        pass


class SequenceAnnotation(models.Model):
    """
    Identify features and sub-components on a DNA or protein sequence.
    """
    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=True,
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
    
    ## follow generic relations doc to create link to any type of component
    ## https://docs.djangoproject.com/en/1.3/ref/contrib/contenttypes/
    ## http://stackoverflow.com/questions/933092/generic-many-to-many-relationships
    ## subComponent = ...
    

class Collection(models.Model):
    """
    Collection of parts
    """
    #: required ID
    displayId = models.CharField('ID', max_length=20,
                    help_text='Identification; Unique within context.', 
                    unique=True)
    
    #: optional short name
    name = models.CharField('Name', max_length=50, blank=True, null=True,
                            help_text='Descriptive name.')
    
    #: required short description -- will be added as first line in SBOL
    #: exports (where this field doesn't exist)
    shortDescription = models.CharField( 'short Description', max_length=200,
        help_text='Very short description for listings')

    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=True,
                help_text='Specify external URI or leave blank for local URI')
    
    #: optional long description
    description = models.TextField( 'Detailed description', blank=True,
        help_text='Use restructured text markup for links, lists and headlines.')
    
    dnaComponents = models.ManyToManyField( DnaComponent, null=True, blank=True,
                                            verbose_name='DNA parts' )

    proteinComponents = models.ManyToManyField( ProteinComponent, 
                                                null=True, blank=True,
                                                verbose_name='Protein parts' )

    chassis = models.ManyToManyField( Chassis, 
                                      null=True, blank=True,
                                      verbose_name="Cell types / Chassis")
    
    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s %s' % (self.displayId, name)
     