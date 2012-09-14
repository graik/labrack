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

from django.contrib.auth.models import User
from django.contrib.auth.models import Group
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
from django.core.files.storage import FileSystemStorage
from django.db import models
from django.db.models import Q
from django.utils.safestring import mark_safe
from django.utils.safestring import SafeUnicode

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
    
    def __unicode__(self):
        return self.displayId + " (" + self.temperature.__str__() \
               + unichr(176) + "C)  [Room: " + self.room + "]"

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'location/%i/' % self.id




class PermissionModel( models.Model ):
    #: Permissions
    created_by = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='%(class)s_created_by')

    owners = models.ManyToManyField(User, null=True, blank=True, 
                                    related_name='%(class)s_owners')
    group_read = models.ManyToManyField(Group, null=True, blank=True, 
                                        related_name='%(class)s_groups_read')
    group_write = models.ManyToManyField(Group, null=True, blank=True, 
                                         related_name='%(class)s_groups_write')
 
   
    def writePermission(self, currentUser):
        
        # Superusers can modify
        if currentUser.is_superuser:
            return True
        
        # Owners can modify
        for user in self.owners.all():
            if user == currentUser:
                return True
            
        # Groups with write can modify
        for group_obj in self.group_write.all():
            for group_member in currentUser.groups.all():
                if group_obj == group_member:
                    return True        
        
        return False
   
    class Meta:
        abstract = True



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
    
    location = models.ForeignKey(Location, related_name='containers')



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




class SampleCollection(PermissionModel):
    """
    SampleCollection to group sample together
    """
    
    name = models.CharField(max_length=25, unique=True)

    description = models.TextField( 'Detailed description', blank=True)
    
    creation_date = models.DateTimeField(auto_now_add=True)
    
    modification_date = models.DateTimeField(auto_now=True)
    

    def __unicode__(self):
        return self.name




class Sample( PermissionModel ):
    """
    Sample describes a single tube or well holding DNA, cells or protein.
    
    The content of the sample (objects derived from Component) is linked in via
    a separate table (SampleContent) together with amount, amountUnit,
    concentration and concentrationUnit. The content objects can be accessed
    via a related manager::
    
        content = sample.sampleContent.objects.all()
        if content[0].content_type.model == u'dnacomponent':
            content[0].concentration = 100.
            content[0].concentrationUnit = 'ng/ul'
       
    or by type of content::
    
       sample.sampleContent.filter( content_type__model=u'dnacomponent' )
       
    We have created shortcut properties that give direct access to all DNA, 
    protein, peptide, chemical or chassis content objects::
    
       dna_content = sample.dnas
       dna_content
       >>>[<DnaComponent: sb0002 - TEV protease site>]
       
       pro_content = sample.proteins
       pep_content = sample.peptides
       che_content = sample.chemicals
       cell_content= sample.chassis
       
    Or all objects regardless of which sub-class of Component they belong to::
    
       sample.components
       >>> [<DnaComponent: sb0002 - TEV protease site>, <ChemicalComponent: bufTE - TE buffer>]
       
    These properties are currently read-only and return a list of standard
    model objects of type DnaComponent, ProteinComponent, etc. If the sample 
    does not have any objects of this type, the property will return an empty
    list.
    """
    
    displayId = models.CharField(max_length=20, 
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
    
    description = models.TextField('Description / Comments', blank=True)
    
    sampleCollection = models.ManyToManyField(SampleCollection, null=True, blank=True)
    
    preparation_date = models.DateField(default=datetime.now())
    
    creation_date = models.DateTimeField(auto_now_add=True)
    
    modification_date = models.DateTimeField(auto_now=True)


    
    class Meta:
        unique_together = ('displayId', 'container')
        ordering = ('container', 'displayId')
    
    # custom properties (shortcuts)
    def _contentObjects( self, model=u'component' ):
        """
        Get Component objects linked through SampleContent entries.
        @param model: unicode, model name
        @return [ Component ] or [] if none
        """
        l = self.sampleContent.filter(content_type__model=model)
        if not l:
            return []
        return [ o.content_object for o in l ]

    @property
    def components( self ):
        """
        Get all component objects linked through SampleContent entries.
        @return: [ Component-derived ] or [] if there is no content
        """
        return self.dnas + self.proteins + self.peptides + self.chassis\
            + self.chemicals
    
    @property
    def dnas( self ):
        """
        DnaComponent objects linked through sample content entries.
        @return [ DnaComponent ], list of DnaComponent objects or [] if none
        """
        return self._contentObjects( u'dnacomponent' )
   
    @property
    def proteins( self ):
        """
        ProteinComponent objects linked through sample content entries.
        @return [ ProteinComponent ], list of ProteinComponent objects or [] if none
        """
        return self._contentObjects(u'proteincomponent')

    @property
    def chassis( self ):
        """
        Chassis objects linked through sample content entries.
        @return [ Chassis ], list of Chassis objects or [] if none
        """
        return self._contentObjects(u'chassis')
    
    @property
    def peptides( self ):
        """
        PeptideComponent objects linked through sample content entries.
        @return [ PeptideComponent ], list of PeptideComponent objects or [] if none
        """
        return self._contentObjects(u'peptidecomponent')
    
    @property
    def chemicals( self ):
        """
        ChemicalComponent objects linked through sample content entries.
        @return [ ChemicalComponent ], list of ChemicalComponent objects or [] if none
        """
        return self._contentObjects(u'chemicalcomponent')
    
    @property
    def sampleType( self ):
        """
        @return str: Either of 'cells', 'DNA', 'protein', 'peptide', 'chemical'
        """
        if self.chassis:
            return u'cells'
        if self.dnas:
            return u'DNA'
        if self.proteins:
            return u'protein'
        if self.peptides:
            return u'peptide'
        if self.chemicals:
            return u'chemical'
        return 'unknown'

    # override Django object methods
    def __unicode__(self):
        return u'%s' % (self.displayId)

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'sample/%i/' % self.id
    
    
    # custom display methods for web interface        
    def qr_code(self):
        """?"""
        data = str(self.displayId + '\n' + self.name + '\n' \
                   + self.preparation_date.__str__() + '\n')
        
        for x in self.samplecontent.all():
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
        
    def strFullContent(self):
        """
        @return str, all content items including amount and concentration
        """
        contentString = u''.encode( 'utf-8' )
        
        for sc in self.sampleContent.all():
            o = sc.content_object
            
            if contentString:
                contentString += u'\n'     ## add line break to previous line

            contentString += unicode(o)

            amount = u'%5.1f %s' %(sc.amount, sc.amountUnit) if sc.amount else ''
            conc =  u'%5.1f %s' % (sc.concentration, sc.concentrationUnit) \
                if sc.concentration else ''
            
            if amount or conc:
                contentString += u' [%s %s]' % (amount, conc) 

        return contentString
        

    def showFullContent(self):
        """
        Show all content items including amount and concentration with URL 
        links.
        """
        contentString = u''.encode( 'utf-8' )
        
        for sc in self.sampleContent.all():
            o = sc.content_object
            
            if contentString:
                contentString += u'<br>'     ## add line break to previous line

            contentString += u'<a href="%s/%s">%s</a>' % (S.admin_root, 
                                                         o.get_relative_url(),
                                                         str(o) )

            amount = u'%5.1f %s' %(sc.amount, sc.amountUnit) if sc.amount else ''
            conc =  u'%5.1f %s' % (sc.concentration, sc.concentrationUnit) \
                if sc.concentration else ''
            
            if amount or conc:
                contentString += u' [%s %s]' % (amount, conc) 

        return contentString

    showFullContent.allow_tags = True  ## don't HTML-escape this string
    showFullContent.short_description = u'Content'
    
    
    def showObject(self, olist, description='', brief=False ):
        """pick first Component object from list and wrap it into URL"""
        o = olist[0]
        item = o.displayId if brief else str(o)
        
        s = u'<a href="%s/%s" title="%s">%s</a>' % \
            (S.admin_root, o.get_relative_url(), description, item )
        
        if len(olist) > 1:
            s += u' + more'
        
        return s
        
        
    def showMainContent(self):
        """
        Only show primary content objects and full content as mouse-over.
        Scenarios:
        Ignore chemicals if there is dna, protein, peptide or chassis;
        
        """
        description = self.strFullContent()
        s = u''
        if self.dnas:
            s = self.showObject( self.dnas, description)
        elif self.proteins:
            s = self.showObject( self.proteins, description )
        elif self.peptides:
            s = self.showObject( self.peptides, description )
            
        if s and self.chassis:
            s += u' in %s' % self.showObject( self.chassis, brief=True ) 
        
        if s:
            return s
        
        if self.chassis:
            return self.showObject( self.chassis, description )
        if self.chemicals:
            return self.showObject( self.chemicals, description )

        return u''
    showMainContent.allow_tags = True  ## don't HTML-escape this string
    showMainContent.short_description = u'Content'
    
    
    
    def strProvenance(self):
        pedigreeString = ""
        for sp in self.samplepedigree.all():
            pedigreeString += sp.sample_source.__str__() + ": "
            if(str(sp.amount) != "None"):
                pedigreeString += ' a[' + str(sp.amount) + ' ' + str(sp.amount_unit) + ']' 
            if(str(sp.concentration) != "None"):
                pedigreeString += ' c[' + str(sp.concentration) + ' ' + str(sp.concentration_unit)  + ']'
            pedigreeString += '\n'
        return pedigreeString
    
    def sampleLinkStr(self):
        linkString = ""
        for sl in self.samplelink.all():
            linkString += sl.linkType + " : " + sl.link + " [" + sl.name + "]"
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
    def showDna(self):
        """filter '(None)' display in admin table"""
        if self.dnas:
            return self.dnas[0]
        return u''
    showDna.short_description = 'DNA'
    showDna.admin_order_field = 'dna'

#    def show_vector(self):
#        """filter '(None)' display in admin table"""
#        if self.vector:
#            return self.vector
#        return u''
#    show_vector.short_description = 'in Vector'
#    show_vector.admin_order_field = 'vector'
#
    def showCell(self):
        """filter '(None)' display in admin table"""
        if self.chassis:
            return self.chassis[0]
        return u''
    showCell.short_description = 'in Cell'
    showCell.admin_order_field = 'cell'
    
    def showCellLink(self):
        """filter '(None)' display in admin table and give link"""
        if self.chassis:
            return self.showObject( self.chassis )
        return u''

    def showProtein(self):
        """filter '(None)' display in admin table"""
        if self.proteins:
            return self.proteins[0]
        return u''
    showProtein.short_description = 'Protein'
    showProtein.admin_order_field = 'protein'
    
    def showSampleType( self ):
        """
        @return: str; either of: 'cells', 'DNA', 'protein', 'peptide', 'chemical'
        """
        return self.sampleType
    showSampleType.short_description = 'Type'

    def showId( self ):
        """
        @return: str; full ID composed of container-sample.
        """
        return u'%s - %s' % (self.container.displayId, self.displayId)
    showId.short_description = 'full ID'
    showId.allow_tags = True  ## don't HTML-escape this string
    
    def showComment( self ):
        """
        @return: str; truncated comment
        """
        if not self.description:
            return u''
        if len(self.description) < 40:
            return unicode(self.description)
        return unicode(self.description[:38] + '..')
    showComment.short_description = 'comments'
#    
##    def get_sequence( self, recenter=0 ):
##        return self.dna.get_sequence( recenter=recenter )
##
##    def get_pretty_sequence( self, recenter=0 ):
##        return self.dna.get_pretty_sequence( recenter=recenter )







class SampleContent(models.Model):
    """
    Helper class to define the content that the sample is made of.
    """
    COMPONENT_LIMITS = {'model__in':('chemicalcomponent', 
                                     'dnacomponent', 
                                     'peptidecomponent', 
                                     'proteincomponent',
                                     'chassis'
                                     )
                        }   

    sample = models.ForeignKey(Sample, related_name='sampleContent')

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





class SampleProvenance(models.Model):
    """
    Define the components that the sample is made of.
    """
    sample = models.ForeignKey(Sample, related_name='sampleProvenance')
    
    sample_source = models.ForeignKey(Sample, null=True, related_name='sampleSource')
    
    PROVENANCE_TYPE = (('parent','parent'), )
    
    provenanceType = models.CharField(max_length=25, choices=PROVENANCE_TYPE)
    
    shortDescription = models.CharField( max_length=50,
                                         help_text='Brief description for tables\
                                         and listings', null=True, blank=True)
    




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





class Component(PermissionModel):
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

 
 
    def formatedUrl(self):
        name = self.name if self.name else ''
        return SafeUnicode("<a href='/admin/labrack/" +self.get_relative_url() \
                           + "'>" + self.displayId + ' - ' +  name + "</a>")

    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s - %s' % (self.displayId, name)
    
   
    class Meta:
        abstract = False
 

    def related_samples( self ):
        """
        """
        r = SampleContent.objects.filter( object_id=self.id, 
                        content_type__model=self.__class__.__name__.lower() )
        s = [content.sample for content in r]
        
        return s
        
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





class PeptideComponent(Component):
    """
    Description of a peptide 'part'.
    """
    sequence = models.TextField( help_text='amino acid sequence', 
                                 blank=True, null=True )
    
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'peptidecomponent/%i/' % self.id

    def size(self):
        """@return int; size in aminoacids"""
        if self.sequence:
            return len( self.sequence )
        return 0
    
    
    
#    def related_samples( self ):
#        """
#        """
#        return self.protein_samples.all()





class ChemicalComponent(Component):
    """
    Description of a chemical.
    """
    #: To define
    
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
    
    #: required short description -- will be added as first line of description
    #: in SBOL exports (where this field doesn't exist)
    shortDescription = models.CharField( 'short Description', max_length=200,
        help_text='Very short description for listings')

    #: uri should be constructed from #displayId if left empty
    uri = models.URLField( blank=True, null=True, unique=False,
                help_text='Specify external URI or leave blank for local URI')
    
    #: optional long description
    description = models.TextField( 'Detailed description', blank=True,
        help_text='Use restructured text markup for links, lists and headlines.')
    
    components = models.ManyToManyField( Component, null=True, blank=True,
                                         verbose_name='parts' )

    
    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s %s' % (self.displayId, name)
