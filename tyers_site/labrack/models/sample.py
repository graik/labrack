from django.db import models
from permissionModel import PermissionModel
from labrack.models.generalmodels import Container
from labrack.models.component import DnaComponent
from labrack.models.component import Chassis
from labrack.models.generalmodels import ContentType
from django.contrib.contenttypes import generic
from labrack.models.unit import Unit


from csvImporter.model import CsvModel
from csvImporter import fields

from django.db.models import Q


from datetime import datetime


## importing custom models

import os


import urllib
import tyers_site.settings as S






class SampleContent(models.Model):
    """
    Helper class to define the content that the sample is made of.
    """
    COMPONENT_LIMITS = {'model__in':('component',
                                     'chemicalcomponent', 
                                     'dnacomponent', 
                                     'peptidecomponent', 
                                     'proteincomponent',
                                     'chassis'
                                     )
                        }   

    sample = models.ForeignKey('Sample', related_name='sampleContent')

    content_type = models.ForeignKey(ContentType, 
                                     verbose_name='Type of content',
                                     limit_choices_to=COMPONENT_LIMITS)

    object_id = models.PositiveIntegerField('Content')

    content_object = generic.GenericForeignKey('content_type', 'object_id')

    solvent = models.CharField('Solvent/Buffer', max_length=100, 
                               blank=True)

    concentration = models.FloatField('Concentration', null=True, blank=True)

    concentrationUnit = models.ForeignKey(Unit, 
                                          verbose_name='Concentration Unit',
                                          related_name='concentrationUnit', 
                                          null=True, blank=True, 
                                          limit_choices_to = {'unitType': 'concentration'})

    amount = models.FloatField('Amount', null=True, blank=True)

    amountUnit = models.ForeignKey(Unit, verbose_name='Amount Unit',
                                   related_name='amountUnit', 
                                   null=True, blank=True, 
                                   limit_choices_to=\
                                   Q(unitType='volume') | Q(unitType='mass') | \
                                   Q(unitType='number')
                                   )

    class Meta:
        app_label = 'labrack'     




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
    class Meta:
        app_label = 'labrack'   


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
    reference_status = models.BooleanField('Reference sample',default=False,
                                           help_text='mark sample as master or reference')


    description = models.TextField('Description / Comments', blank=True)


    sampleCollection = models.ManyToManyField(SampleCollection, null=True, blank=True)

    preparation_date = models.DateField(default=datetime.now(), verbose_name="date")

    creation_date = models.DateTimeField(auto_now_add=True)

    modification_date = models.DateTimeField(auto_now=True)


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
    
    derivedFrom = models.ForeignKey('self', verbose_name='Derived From', blank=True, null=True, related_name='derivedFromSample')

    TYPE_CHOICES = (('parent', 'parent'), 
                    ('mixture', 'mixture'),
                    )
    provenanceType = models.CharField( max_length=30, choices=TYPE_CHOICES, null=True, blank=True, 
                                       default='') 
    historyDescription = models.CharField('Comment', null=True,blank=True,max_length = 255)    

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
        return 'dnasample/%i/' % self.id


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
        try:      
            o = olist[0]
            item = o.displayId if brief else str(o)

            s = u'<a href="%s/%s" title="%s">%s</a>' % \
                (S.admin_root, o.get_relative_url(), description, item )

            if len(olist) > 1:
                s += u' + more'        
            return s
        except:
            s = "---"
            print 'Unexpected error'



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

    def showPeptide(self):
        """filter '(None)' display in admin table"""
        if self.peptides:
            return self.peptides[0]
        return u''
    showPeptide.short_description = 'Peptide'
    showPeptide.admin_order_field = 'peptide'


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

    def is_reference(self):
        return bool(self.reference_status)

    showComment.short_description = 'comments'

    class Meta:    
        app_label = 'labrack' 
        unique_together = ('displayId', 'container')
        ordering = ('container', 'displayId')    
##    def get_sequence( self, recenter=0 ):
##        return self.dna.get_sequence( recenter=recenter )
##
##    def get_pretty_sequence( self, recenter=0 ):
##        return self.dna.get_pretty_sequence( recenter=recenter )

class CSVSample(CsvModel):
    DisplayId = fields.CharField()
    name = fields.CharField()
    container = fields.CharField()
    status = fields.IntegerField()

    PreparationDate = fields.CharField()
    NumberOfAliquots = fields.CharField()
    isReference = fields.CharField()
    Description	= fields.CharField()
    SampleCollection = fields.CharField()

    Chemical_DisplayID = fields.CharField()
    Chemical_SolventBuffer = fields.CharField()
    Chemical_Concentration = fields.CharField()
    Chemical_Concentration_Unit = fields.CharField()
    Chemical_Amount = fields.CharField()
    Chemical_Amount_Unit = fields.CharField()

    DNA_DisplayID = fields.CharField()
    DNA_SolventBuffer = fields.CharField()
    DNA_Concentration = fields.CharField()
    DNA_Concentration_Unit = fields.CharField()
    DNA_Amount = fields.CharField()
    DNA_Amount_Unit = fields.CharField()

    Peptide_DisplayID = fields.CharField()
    Peptide_SolventBuffer = fields.CharField()
    Peptide_Concentration = fields.CharField()
    Peptide_Concentration_Unit = fields.CharField()
    Peptide_Amount = fields.CharField()
    Peptide_Amount_Unit = fields.CharField()

    Protein_DisplayID = fields.CharField()
    Protein_SolventBuffer = fields.CharField()
    Protein_Concentration = fields.CharField()
    Protein_Concentration_Unit = fields.CharField()
    Protein_Amount = fields.CharField()
    Protein_Amount_Unit = fields.CharField()

    class Meta: 
        delimiter = ";"
        has_header = True


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


    class Meta:    
        app_label = 'labrack' 


class DnaSample(Sample):
    #, help ='Either select an existing DNA Construct or fill the dna description to create a new one'
    dnaConstruct = models.ForeignKey(DnaComponent,verbose_name='DNA Construct', blank=True, null=True, related_name='dnaSample',help_text='Either select an existing DNA Construct or fill the dna description to create a new one')
    
    
    inChassis = models.ForeignKey(Chassis, verbose_name='in Cell', blank=True, null=True, related_name='chassisSample')

    
    class Meta:
        app_label = 'labrack'
        verbose_name = 'DNA Sample'
        abstract = False
        
    def related_samples( self ):
        """
        """
        
        r = DnaSample.objects.filter( dnaConstruct = self.dnaConstruct.id ) 
     
        return r
    
    def get_relative_url(self):
            """
            Define standard relative URL for object access in templates
            """
            return 'dnasample/%i/' % self.id   
        
    def Content(self):
        return self.dnaConstruct.formatedUrl();    
    
    def inCell(self):
            return self.chassis.formatedUrl();   
    
    
    
class ChassisSample(Sample):
    #, help ='Either select an existing DNA Construct or fill the dna description to create a new one'
    chassis = models.ForeignKey(Chassis,verbose_name='Chassis', null=True, blank=True, related_name='Cell')
        
    class Meta:
        app_label = 'labrack'
        verbose_name = 'Cell Sample'
        abstract = False
        
    def related_samples( self ):
        """
        """
        r = ChassisSample.objects.filter( chassis = self.chassis.id ) 
     
        return r   
    
    def get_relative_url(self):
            """
            Define standard relative URL for object access in templates
            """
            return 'chassissample/%i/' % self.id      