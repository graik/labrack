from django.db import models
from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from datetime import datetime
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic


# Create your models here.



################################################################################################################
class Location(models.Model):
    """
    A location (fridge, freezer, shelf) where container are stored
    """
    displayId = models.CharField('Location (ID)', max_length=20, unique=True, help_text='Unique identifier')
    shortDescription = models.CharField('Short description', max_length=200, blank=True, help_text='Brief description of fridge, freezer or shelf')
    temperature = models.FloatField('Temperature', help_text= unichr(176) + 'C')
    room = models.CharField('Room', max_length=20)
    creation_date = models.DateTimeField(auto_now_add=True)
    modification_date = models.DateTimeField(auto_now=True)
    def __unicode__(self):
        return self.displayId + " (" + self.temperature.__str__() + unichr(176) + "C)  [Room: " + self.room + "]"

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'location/%i/' % self.id



################################################################################################################
class Container( models.Model ):
    """
    A container holding several physical samples of nucleic acids, proteins or other stuff.
    """

    STORAGE_CONTAINER_TYPES= (
        ('96-well-plate', '96 well plate'),
        ('384-well-plate','384 well plate'),
        ('box', 'Freezer box'),
        ('other', 'other' ) )
    
    displayId = models.CharField('Container (ID)', max_length=20, unique=True, help_text='Unique identifier. E.g. D001 or C012 (D...DNA, C...Cells)')
    shortDescription = models.CharField('Short description', max_length=200, blank=True, help_text='Brief description for tables and listings')
    containerType = models.CharField('Type of container', max_length=30, choices=STORAGE_CONTAINER_TYPES )
    location = models.ForeignKey(Location, related_name='containers')

    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True, related_name='container_created_by')
    owners = models.ManyToManyField(User, null=True, blank=True, related_name='container_owners')
    group_read = models.ManyToManyField(Group, null=True, blank=True, related_name='container_groups_read')
    group_write = models.ManyToManyField(Group, null=True, blank=True, related_name='container_groups_write')

    #: optional long description
    description = models.TextField( 'Detailed description', blank=True)
    creation_date = models.DateTimeField(auto_now_add=True)
    modification_date = models.DateTimeField(auto_now=True)



    def __unicode__( self ):
        return "%s" % self.displayId
    
    class Meta:
        ordering = ('displayId',)

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'container/%i/' % self.id
    

#    def next_well( self ):
#        """
#        Try to guess next free well.
#        """
#        r = ''
#
#        if self.container_type == 'box':
#            try:
#                wells = [ int(s.well) for s in self.samples.all() ]
#                wells.sort()
#                r = str( wells[-1] + 1 )
#            except:
#                r = "%02i" % self.samples.count() + 1
#
#        if self.container_type in ['96-well-plate', '384-well-plate']:
#            try:
#                wells = [ s.well for s in self.samples.all() ]
#                wells.sort()
#                row = wells[-1][0]
#                col = wells[-1][1:] + 1
#                r = "%s%02i" % (row, col)
#            except:
#                r = "%02i" % self.samples.count() + 1
#
#        return r
    

    




################################################################################################################    
class Unit(models.Model):
    """
    Unit for amount, concentration, volume
    """
    name = models.CharField(max_length=10)
    type = models.CharField(max_length=25)
    
    def __unicode__(self):
        return self.name




################################################################################################################
class Sample( models.Model ):
    """
    Sample describes a single tube or well holding DNA, cells or protein.
    """

    STORAGE_TYPES      = ( ('dna', 'DNA'), ('cells','cells'), ('protein','protein') )
    CONCENTRATION_UNITS= ( ('mg/l', 'mg / l'), ('mol/l','mol / l'), 
                           ('umol/l',u'\u00B5' + u'mol / l') )

    
    displayId = models.CharField('Sample (ID)', max_length=20, help_text='Label or well position. Must be unique within container.' )
    shortDescription = models.CharField('Short description', max_length=200, blank=True, help_text='Brief description for tables and listings')
    container = models.ForeignKey( Container, related_name='samples' )  #: link to a single container
    aliquotnb = models.PositiveIntegerField('Aliquot number', null=True, blank=True)
    description = models.TextField( 'Detailed description', blank=True)
    preparation_date = models.DateTimeField(default=datetime.now())
    creation_date = models.DateTimeField(auto_now_add=True)
    modification_date = models.DateTimeField(auto_now=True)


    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True, related_name='sample_created_by')
    owners = models.ManyToManyField(User, null=True, blank=True, related_name='sample_owners')
    group_read = models.ManyToManyField(Group, null=True, blank=True, related_name='sample_groups_read')
    group_write = models.ManyToManyField(Group, null=True, blank=True, related_name='sample_groups_write')

    
    class Meta:
        unique_together = ('displayId', 'container')
        ordering = ('container', 'displayId')
    

    def __unicode__( self ):
        return self.container.displayId + ' / ' + self.displayId

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'sample/%i/' % self.id
    
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

#    def related_samples( self ):
#        """
#        @return: QuerySet of samples with same DNA content
#        """
#        r = []
#
#        if self.dna:
#            r = Sample.objects.filter( dna=self.dna )
#        elif self.vector:
#            ## only consider vector if there is no insert:
#            r = Sample.objects.filter( vector=self.vector )
#
#        elif self.protein:
#            r = Sample.objects.filter( protein=self.protein )
#
#        elif self.cell:
#            r = Sample.objects.filter( cell=self.cell )
#            
#        if r:
#            return r.exclude( pk=self.pk )
#        return r
#    
#
#    def show_dna(self):
#        """filter '(None)' display in admin table"""
#        if self.dna:
#            return self.dna
#        return u''
#    show_dna.short_description = 'DNA'
#    show_dna.admin_order_field = 'dna'
#
#    def show_vector(self):
#        """filter '(None)' display in admin table"""
#        if self.vector:
#            return self.vector
#        return u''
#    show_vector.short_description = 'in Vector'
#    show_vector.admin_order_field = 'vector'
#
#    def show_cell(self):
#        """filter '(None)' display in admin table"""
#        if self.cell:
#            return self.cell
#        return u''
#    show_cell.short_description = 'in Cell'
#    show_cell.admin_order_field = 'cell'
#    
#    def show_protein(self):
#        """filter '(None)' display in admin table"""
#        if self.protein:
#            return self.protein
#        return u''
#    show_protein.short_description = 'Protein'
#    show_protein.admin_order_field = 'protein'
#    
#    def show_sampleType( self ):
#        """
#        @return: str; either of: , 'cells', 'DNA', 'protein' or 'unknown'
#        """
#        if self.cell:
#            return 'cells'
#        if self.dna or self.vector:
#            return 'DNA'
#        if self.protein:
#            return 'protein'
#        return 'unknown'
#    show_sampleType.short_description = 'Type'
#
#    def show_Id( self ):
#        """
#        @return: str; full ID composed of container-sample.
#        """
#        return u'%s - %s' % (self.container.displayId, self.displayId)
#    show_Id.short_description = 'full ID'
#    show_Id.allow_tags = True  ## don't HTML-escape this string
#    
#    def show_concentration( self ):
#        """
#        @return: str; concentration + unit
#        """
#        return u'%3.2f %s' % (self.concentration, self.concentrationUnit)
#    show_concentration.short_description = 'Concentration'
#    show_concentration.admin_order_field = 'concentration'
#    
#    def show_comments( self ):
#        """
#        @return: str; truncated comment
#        """
#        if not self.comments:
#            return u''
#        if len(self.comments) < 40:
#            return unicode(self.comments)
#        return unicode(self.comments[:38] + '..')
#    show_comments.short_description = 'comments'
#    
##    def get_sequence( self, recenter=0 ):
##        return self.dna.get_sequence( recenter=recenter )
##
##    def get_pretty_sequence( self, recenter=0 ):
##        return self.dna.get_pretty_sequence( recenter=recenter )

    


################################################################################################################

relation_limits = {'model__in':('dnacomponent', 'proteincomponent')}
class SampleComponent(models.Model):
    sample = models.ForeignKey(Sample)
    content_type = models.ForeignKey(ContentType, limit_choices_to=relation_limits)
    object_id = models.PositiveIntegerField()
    content_object = generic.GenericForeignKey('content_type', 'object_id')
    amount = models.FloatField('Amount/Volume', null=True, blank=True)
    amount_unit = models.ForeignKey(Unit, related_name='amount_unit', null=True, blank=True)
    concentration = models.FloatField('Concentration', null=True, blank=True)
    concentration_unit = models.ForeignKey(Unit, related_name='concentration_unit', null=True, blank=True)


################################################################################################################
class ComponentType( models.Model ):
    """
    Helper class for classifying parts.
    Following SBOL, each type should correspond to a Sequence Ontology term.
    """
    
    #: required
    uri = models.URLField( unique=False,
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






################################################################################################################
class Component(models.Model):
    """
    Abstract base class for DNA and protein 'parts' as well as host cells.
    """
    
    STATUS_CHOICES = ( ('available', 'available'),
                       ('planning', 'planning'),
                       ('under_construction', 'under construction'),
                       ('abandoned', 'abandoned'))
    
    #: required ID
    displayId = models.CharField('ID', max_length=20,
                    help_text='Identification; Unique within collection.', 
                    unique=True)
    
    #: optional short name
    name = models.CharField('Name', max_length=50, blank=True, null=True,
                            help_text='descriptive name (e.g. "TEV protease")')
    
    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=False,
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
    
    #: non-SBOL compliant status
    status = models.CharField( max_length=30,
                               choices=STATUS_CHOICES,
                               default='planning')
    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s - %s' % (self.displayId, name)
    
   
    class Meta:
        abstract = True
 
    def isavailable(self):
        return self.samples.count() > 0
    isavailable.short_description = 'available'
    isavailable.boolean = True
    
    def related_samples( self ):
        """
        """
        from django.db.models import Q
        q = Q(dna=self) | Q(cell=self) | Q(vector=self) | Q(protein=self)
        
        r = Sample.objects.filter( q )
        return r
    
    



################################################################################################################
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
    
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'dnacomponent/%i/' % self.id

    def size(self):
        """@return int; size in nucleotides"""
        if self.sequence:
            return len( self.sequence )
        return 0

    def isavailable_dna(self):
        """True if there is a DNA-only sample containing this DnaComponent"""
        return self.dna_samples.count() > 0
    isavailable_dna.short_description = 'DNA'
    isavailable_dna.boolean = True
    
    def isavailable_cells(self):
        """
        True if there is a cell stock (DNA + cells) sample for this 
        DnaComponent
        """
        samples = self.dna_samples.all()

        for s in samples:
            if s.cell:
                return True
        return False
    isavailable_cells.short_description = 'Cells'
    isavailable_cells.boolean = True
    
    def related_samples( self ):
        """
        """
        return self.dna_samples.all()

    def show_translatesTo(self):
        """filter '(None)' display in admin table"""
        if self.translatesTo:
            return self.translatesTo
        return u''
    show_translatesTo.short_description = 'translates to'
    show_translatesTo.admin_order_field = 'translatesTo'
     
    def show_optimizedFor(self):
        """filter '(None)' display in admin table"""
        if self.optimizedFor:
            return self.optimizedFor
        return u''
    show_optimizedFor.short_description = 'optimized for'
    show_optimizedFor.admin_order_field = 'optimizedFor'
     

    
################################################################################################################
class SelectiveMarker( models.Model ):
    """
    Describes an Antibiotic or similar resistence marker
    """
    #: optional short name
    name = models.CharField('Name', max_length=50, blank=True, null=True,
                            help_text='descriptive name (e.g. "TEV protease")')
    
    #: required short description -- will be added as first line in SBOL
    #: exports (where this field doesn't exist)
    shortDescription = models.CharField( 'short Description', max_length=200,
        help_text='very short description for listings')

    #: optional long description
    description = models.TextField( 'Detailed description', 
                                    blank=True, null=True )

    def __unicode__(self):
        return self.name

##    def get_absolute_url(self):
##        """Define standard URL for object.get_absolute_url access in templates """
##        return APP_URL+'/selectivemarker/%i/' % self.id

    class Meta:
        ordering = ('name',)
    





################################################################################################################    
class VectorDnaComponent(DnaComponent):
    """
    Description of vector backbone. So far identical to DnaComponent.
    """

    #: [optional] link to resistence or similar marker
    marker = models.ManyToManyField( SelectiveMarker,
                                     verbose_name='selective marker(s)',
                                     null=True, blank=True)

    class Meta(Component.Meta):
        pass
    
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'vectordnacomponent/%i/' % self.id

    def related_samples( self ):
        """
        """
        return self.vector_samples.all()
    




################################################################################################################    
class ProteinComponent(Component):
    """
    Description of a amino-acid encoded protein 'part'.
    """
    #: optional sequence
    sequence = models.TextField( help_text='amino acid sequence', 
                                 blank=True, null=True )

    class Meta(Component.Meta):
        pass

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'proteincomponent/%i/' % self.id

    def size(self):
        """@return int; size in aminoacids"""
        if self.sequence:
            return len( self.sequence )
        return 0

    def related_samples( self ):
        """
        """
        return self.protein_samples.all()




################################################################################################################
class Chassis(Component):
    """
    Description of a host system. Usually this will be a cell type or bacterial
    strain.
    """

    class Meta(Component.Meta):
        pass

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'chassis/%i/' % self.id

    def related_samples( self ):
        """
        """
        return self.cell_samples.all()




################################################################################################################
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
    
    ## follow generic relations doc to create link to any type of component
    ## https://docs.djangoproject.com/en/1.3/ref/contrib/contenttypes/
    ## http://stackoverflow.com/questions/933092/generic-many-to-many-relationships
    ## subComponent = ...
    




################################################################################################################
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
    uri = models.URLField( blank=True, null=True, unique=False,
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
     
     
     
     
################################################################################################################

