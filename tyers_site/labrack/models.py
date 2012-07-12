from django.db import models
from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from datetime import datetime
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
import urllib
import tyers_site.settings as settings

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
    
    def __unicode__(self):
        return self.displayId + " (" + self.temperature.__str__() \
               + unichr(176) + "C)  [Room: " + self.room + "]"

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'location/%i/' % self.id



class Container( models.Model ):
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
    
    location = models.ForeignKey(Location, related_name='containers')


    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='container_created_by')
    owners = models.ManyToManyField(User, null=True, blank=True, 
                                    related_name='container_owners')
    group_read = models.ManyToManyField(Group, null=True, blank=True, 
                                        related_name='container_groups_read')
    group_write = models.ManyToManyField(Group, null=True, blank=True, 
                                         related_name='container_groups_write')

    #: optional long description
    description = models.TextField( 'Detailed description', blank=True)
    
    creation_date = models.DateTimeField('Created at', auto_now_add=True)
    
    modification_date = models.DateTimeField('Modified at', auto_now=True)


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

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'container/%i/' % self.id
    

class Unit(models.Model):
    """
    Unit for amount, concentration, volume
    """
    
    UNIT_TYPE = (('mass','mass'), 
                 ('volume','volume'), 
                 ('concentration','concentration'),
                 ('number','number'),
                 ('other','other'))
    
    name = models.CharField(max_length=10)
    
    conversion = models.FloatField('Conversion Factor', blank=True, null=True,
                            help_text='Factor for conversion to SI unit')
    
    unitType = models.CharField(max_length=25, choices=UNIT_TYPE)
    
    
    
    def __unicode__(self):
        return self.name


class Sample( models.Model ):
    """
    Sample describes a single tube or well holding DNA, cells or protein.
    """
    
    displayId = models.CharField(max_length=20, 
        help_text='Label or well position. Must be unique within container.')

    name = models.CharField('Name', max_length=200, blank=True, 
                        help_text='Informative name for tables and listings')
    
    #: link to a single container
    container = models.ForeignKey( Container, related_name='samples' )

    aliquotNr = models.PositiveIntegerField('Number of aliquots', 
                                            null=True, blank=True)

    empty = models.BooleanField('Empty', default=False)
    
    description = models.TextField('Description / Comments', blank=True)
    
    #: date AND time is perhaps a bit overdoing it --Raik
    preparation_date = models.DateTimeField(default=datetime.now())
    
    creation_date = models.DateTimeField(auto_now_add=True)
    
    modification_date = models.DateTimeField(auto_now=True)


    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='sample_created_by')
    owners = models.ManyToManyField(User, null=True, blank=True, 
                                    related_name='sample_owners')
    group_read = models.ManyToManyField(Group, null=True, blank=True, 
                                        related_name='sample_groups_read')
    group_write = models.ManyToManyField(Group, null=True, blank=True, 
                                         related_name='sample_groups_write')

    
    class Meta:
        unique_together = ('displayId', 'container')
        ordering = ('container', 'displayId')
    

    def __unicode__(self):
        return u'%s - %s' % (self.displayId, self.name)

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'sample/%i/' % self.id
    
        
    def qr_code(self):
        """?"""
        data = str(self.displayId + '\n' + self.name + '\n' \
                   + self.preparation_date.__str__() + '\n')
        
        for x in self.sameplecontent.all():
            data += x.content_object.__str__()
            
            if(str(x.amount) != "None"):
                data += '[' + str(x.amount) + ' ' + str(x.amount_unit) + ']' 

            if(str(x.concentration) != "None"):
                data += '[' + str(x.concentration) + ' ' \
                    + str(x.concentration_unit)  + ']'

            data += '\n'
            
        data += '\n' 
        return urllib.quote(data)
        #return data
        
        
    def sampleContentStr(self):
        contentString = ""
        for sc in self.samplecontent.all():
            contentString += sc.content_object.__str__() + ": "
            if(str(sc.amount) != "None"):
                contentString += ' a[' + str(sc.amount) + ' ' + str(sc.amount_unit) + ']' 
            if(str(sc.concentration) != "None"):
                contentString += ' c[' + str(sc.concentration) + ' ' + str(sc.concentration_unit)  + ']'
            contentString += '\n'
        return contentString
    
    
    def samplePedigreeStr(self):
        pedigreeString = ""
        for sp in self.sameplepedigree.all():
            pedigreeString += sp.sample_source.__str__() + ": "
            if(str(sp.amount) != "None"):
                pedigreeString += ' a[' + str(sp.amount) + ' ' + str(sp.amount_unit) + ']' 
            if(str(sp.concentration) != "None"):
                pedigreeString += ' c[' + str(sp.concentration) + ' ' + str(sp.concentration_unit)  + ']'
            pedigreeString += '\n'
        return pedigreeString
    
    def sampleLinkStr(self):
        linkString = ""
        for sl in self.sameplelink.all():
            linkString += sl.type + " : " + sl.link + " [" + sl.name + "]"
            linkString += '\n'
        return linkString
    
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
    def show_Id( self ):
        """
        @return: str; full ID composed of container-sample.
        """
        return u'%s - %s' % (self.container.displayId, self.displayId)
    show_Id.short_description = 'full ID'
    show_Id.allow_tags = True  ## don't HTML-escape this string
    
    def show_comments( self ):
        """
        @return: str; truncated comment
        """
        if not self.description:
            return u''
        if len(self.description) < 40:
            return unicode(self.description)
        return unicode(self.description[:38] + '..')
    show_comments.short_description = 'comments'
#    
##    def get_sequence( self, recenter=0 ):
##        return self.dna.get_sequence( recenter=recenter )
##
##    def get_pretty_sequence( self, recenter=0 ):
##        return self.dna.get_pretty_sequence( recenter=recenter )

class Project(models.Model):
    """
    Project to group sample together
    """
    
    name = models.CharField(max_length=25, unique=True)

    shortDescription = models.CharField('Short description', max_length=200,
                                        blank=True, help_text='Brief description\
                                        for tables and listings')
    
    description = models.TextField( 'Detailed description', blank=True)
    
    
    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True,
                                   related_name='project_created_by')
    
    owners = models.ManyToManyField(User, null=True, blank=True,
                                    related_name='project_owners')
    
    group_read = models.ManyToManyField(Group, null=True, blank=True,
                                        related_name='project_groups_read')
    
    group_write = models.ManyToManyField(Group, null=True, blank=True,
                                         related_name='project_groups_write')

    creation_date = models.DateTimeField(auto_now_add=True)
    
    modification_date = models.DateTimeField(auto_now=True)
    

    def __unicode__(self):
        return self.name


class SampleLink(models.Model):
    """
    Link to file system or url
    """
    
    sample = models.ForeignKey(Sample, related_name='sameplelink')
    
    LINK_TYPE = (('filesystem','File System'), ('url','URL'))
    
    type = models.CharField(max_length=25, choices=LINK_TYPE)
    
    link = models.CharField(max_length=1000, help_text='Link to file system or url')
    
    shortDescription = models.CharField( max_length=50,
                                         help_text='Brief description for tables\
                                         and listings', null=True, blank=True)
    
    
    
    def __unicode__(self):
        if self.type == 'url':
            from django.utils.safestring import SafeUnicode
            return SafeUnicode("<a href='" + self.link + "'>" + self.link + "</a>")
        return self.link


from django.db.models import Q

class SampleContent(models.Model):
    """
    Helper class to define the content that the sample is made of.
    """
    COMPONENT_LIMITS = {'model__in':('chemicalcomponent', 
                                     'dnacomponent', 
                                     'peptidecomponent', 
                                     'proteincomponent',
    ## I think samples within samples would really complicate our life...--Raik
                                     'sample')}   

    sample = models.ForeignKey(Sample, related_name='samepleContent')

    content_type = models.ForeignKey(ContentType, 
                                     limit_choices_to=COMPONENT_LIMITS)

    object_id = models.PositiveIntegerField()

    content_object = generic.GenericForeignKey('content_type', 'object_id')

    amount = models.FloatField('Amount', null=True, blank=True)

    amountUnit = models.ForeignKey(Unit, related_name='amountUnit', 
                    null=True, blank=True, 
                    limit_choices_to=\
                        Q(unitType='volume') | Q(unitType='mass') | \
                        Q(unitType='number')
                            )

    concentration = models.FloatField('Concentration', null=True, blank=True)

    concentrationUnit = models.ForeignKey(Unit, 
                         related_name='concentrationUnit', 
                         null=True, blank=True, 
                         limit_choices_to = {'unitType': 'concentration'})


class SamplePedigree(models.Model):
    """
    Define the components that the sample is made of.
    """
    sample = models.ForeignKey(Sample, related_name='sameplepedigree')
    
    sample_source = models.ForeignKey(Sample, null=True)
    
    amount = models.FloatField('Amount/Volume', null=True, blank=True)
    
    amount_unit = models.ForeignKey(Unit, related_name='%(class)s_amount_unit', null=True,
                                    blank=True, limit_choices_to = {'type': 'amount'})
    
    concentration = models.FloatField('Concentration', null=True, blank=True)
    
    concentration_unit = models.ForeignKey(Unit, related_name='%(class)s_concentration_unit',
                                           null=True, blank=True,
                                           limit_choices_to = {'type': 'concentration'})




class ComponentType( models.Model ):
    """
    Helper class for classifying parts.
    Following SBOL, each type should in theory correspond to a Sequence Ontology term.
    """
    
    uri = models.URLField( unique=False, blank=True, null=True,
                           help_text='Typically a sequence ontology URI, example: http://purl.obolibrary.org/obo/SO_0000167' )

    name = models.CharField('Name', max_length=200, 
                        help_text='Informative name')
    
    #: required ? directional relationship to parent type or types
    subTypeOf = models.ManyToManyField('self', symmetrical=False, blank=True, 
                                       null=True, related_name='subTypes')

    def __unicode__( self ):
        return unicode(self.name)



class Component(models.Model):
    """
    Abstract base class for cells, nucleic acids, proteins, and chemicals.
    """
    
    STATUS_CHOICES = ( ('available', 'available'),
                       ('planning', 'planning'),
                       ('under_construction', 'under construction'),
                       ('abandoned', 'abandoned'))
    
    displayId = models.CharField('display ID', max_length=20, unique=True, 
                                 help_text='Unique identification')

    name = models.CharField('Name', max_length=200, blank=True, 
                        help_text='Descriptive name (e.g. "TEV protease")')

    description = models.TextField( 'Detailed description', blank=True)
    
    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=False, 
                help_text='Specify external URI or leave blank for local URI') 

    #: non-SBOL compliant status    
    status = models.CharField( max_length=30, choices=STATUS_CHOICES, 
                               default='planning')
    
    componentType = models.ManyToManyField(ComponentType, 
                                            blank=True, null=True, 
                                            verbose_name='Part type', 
                                    help_text='Classification of this part.')

    variantOf = models.ManyToManyField( 'self', symmetrical=False, 
                                        blank=True, null=True, 
                                        verbose_name='Variant of', 
            help_text='Specify part(s) this part is derived from, if any.' )
    
    abstract = models.BooleanField( 'Abstract Part', default=False, 
        help_text='Entry only serves as container to organize related parts.')
    
    annotations = models.ManyToManyField( 'SequenceAnnotation', 
                                          verbose_name='Annotations', 
                                          blank=True, null=True )
    
    creation_date = models.DateTimeField(auto_now_add=True, null=True)
    
    modification_date = models.DateTimeField(auto_now=True, null=True)

 
    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='%(class)s_created_by')
    owners = models.ManyToManyField(User, null=True, blank=True, 
                                    related_name='%(class)s_owners')
    group_read = models.ManyToManyField(Group, null=True, blank=True, 
                                        related_name='%(class)s_groups_read')
    group_write = models.ManyToManyField(Group, null=True, blank=True, 
                                         related_name='%(class)s_groups_write')
 
    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s - %s' % (self.displayId, name)
    
   
    class Meta:
        abstract = True
 
#    def isavailable(self):
#        return self.samples.count() > 0
#    isavailable.short_description = 'available'
#    isavailable.boolean = True
#    
#    def related_samples( self ):
#        """
#        """
#        from django.db.models import Q
#        q = Q(dna=self) | Q(cell=self) | Q(vector=self) | Q(protein=self)
#        
#        r = Sample.objects.filter( q )
#        return r


class DnaComponent(Component):
    """
    Description of a stretch of DNA or RNA.
    """
    #: optional sequence
    sequence = models.TextField( help_text='nucleotide sequence', 
                                 blank=True, null=True )
    
    translatesTo = models.ForeignKey( 'ProteinComponent', blank=True, null=True, 
                                      related_name='encodedBy', 
                        help_text='Protein part this sequence translates to' )

    optimizedFor = models.ForeignKey( 'Chassis', blank=True, null=True )

    class Meta(Component.Meta):
        pass
    
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'nucleicacidcomponent/%i/' % self.id

    def size(self):
        """@return int; size in nucleotides"""
        if self.sequence:
            return len( self.sequence )
        return 0

##    def isavailable_dna(self):
##        """True if there is a DNA-only sample containing this NucleicAcidComponent"""
##        return self.dna_samples.count() > 0
##    isavailable_dna.short_description = 'DNA'
##    isavailable_dna.boolean = True
    
##    def isavailable_cells(self):
##        """
##        True if there is a cell stock (DNA + cells) sample for this 
##        NucleicAcidComponent
##        """
##        samples = self.dna_samples.all()
##
##        for s in samples:
##            if s.cell:
##                return True
##        return False
##    isavailable_cells.short_description = 'Cells'
##    isavailable_cells.boolean = True
    
##    def related_samples( self ):
##        """
##        """
##        return self.dna_samples.all()

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


class ProteinComponent(Component):
    """
    Description of a protein 'part'.
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

#    def related_samples( self ):
#        """
#        """
#        return self.protein_samples.all()


class PeptideComponent(ProteinComponent):
    """
    Description of a peptide 'part'.
    """
    
    class Meta(ProteinComponent.Meta):
        pass

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'peptidecomponent/%i/' % self.id

#    def related_samples( self ):
#        """
#        """
#        return self.protein_samples.all()



class ChemicalComponent(Component):
    """
    Description of a chemical.
    """
    #: To define
    
    class Meta(Component.Meta):
        pass

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'chemicalcomponent/%i/' % self.id


   
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
