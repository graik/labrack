from datetime import datetime
import os

from django.db import models
from django.utils.safestring import mark_safe
from django.utils.safestring import SafeUnicode 
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse

from Bio import SeqIO
from Bio.Seq import Seq

from generalmodels import DnaSequenceAnnotation
from usermixin import UserMixinModel
from tyers_site import settings   

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
    subTypeOf = models.ManyToManyField('self', symmetrical=False, blank=True, 
                                       null=True, related_name='subTypes')
    class Meta:
        app_label = 'labrack' 
        verbose_name = 'DNA PartType'
        abstract = False

class ProteinComponentType( ComponentType ):
    #: required ? directional relationship to parent type or types
    subTypeOf = models.ManyToManyField('self', symmetrical=False, blank=True, 
                                       null=True, related_name='subTypes')
    class Meta:
        app_label = 'labrack'         
        verbose_name = 'Protein PartType'
        abstract = False        

class ChemicalComponentType( ComponentType ):
    #: required ? directional relationship to parent type or types
    subTypeOf = models.ManyToManyField('self', symmetrical=False, blank=True, 
                                       null=True, related_name='subTypes')
    class Meta:
        app_label = 'labrack'   
        verbose_name = 'Chemical PartType'
        abstract = False

class ChassisComponentType( ComponentType ):
    #: required ? directional relationship to parent type or types
    subTypeOf = models.ManyToManyField('self', symmetrical=False, blank=True, 
                                       null=True, related_name='subTypes')
    class Meta:
        app_label = 'labrack'   
        verbose_name = 'Cell PartType'
        abstract = False   

class PeptideComponentType( ComponentType ):
    #: required ? directional relationship to parent type or types
    subTypeOf = models.ManyToManyField('self', symmetrical=False, blank=True, 
                                       null=True, related_name='subTypes')
    class Meta:
        app_label = 'labrack'   
        verbose_name = 'Peptide PartType'
        abstract = False         



        
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

    displayId = models.CharField('ID', max_length=20, unique=True, 
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
    
    def show_component_type(self):
            """
            """
            ret = ""
            comDna= DnaComponent.objects.get(id = self.id)
            try:                       
                for t in comDna.componentType.all():
                    if ret=="":
                        ret += t.name
                    else:
                        ret += ","+t.name
            except Exception, err:
                print err        
            #s = [content.sequenceannotation for content in r]
    
            return ret   
    show_component_type.short_description = 'Type'    

   

    def __unicode__(self):
        name = self.name if self.name else ''
        return u'%s - %s' % (self.displayId, name)


    class Meta:
        abstract = False


    def related_samples( self ):
        """
        """
        import labrack.models as M
        r = M.SampleContent.objects.filter( object_id=self.id, 
                                            content_type__model=self.__class__.__name__.lower() )
        s = [content.sample for content in r]

        return s


    def related_annotations( self ):
        """
        """
        r = DnaSequenceAnnotation.objects.filter( subComponent=self.id).order_by('bioStart')

        return r    

    def get_dna_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'component/%i/' % self.id
    
    
    def related_seq( self ):
        """
        """
        gb_features = ''

        gb_file = settings.MEDIA_ROOT+"/"+os.path.normpath(self.GenBankfile.name)
        for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
            gb_features += gb_record.seq.tostring()
        return gb_features



    def related_file_annotation( self ):

        gb_file = settings.MEDIA_ROOT+"/"+os.path.normpath(self.GenBankfile.name)
        gb_features = ""
        for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
            # now do something with the record
            for ind in xrange(len(gb_record.features)) :
                gb_features += '\n'+ repr(gb_record.features[ind].type) + " Location start : "+ repr(gb_record.features[ind].location._start.position) + " Location end : "+ repr(gb_record.features[ind].location._end.position)
            return gb_features


    def number_related_samples( self ):
        """
        """
        import labrack.models as M


        r = M.DnaSample.objects.filter( dnaConstruct=self.id )
        s = r.count()      

        return s
    number_related_samples.short_description = 'Samples'

    def number_related_annotations( self ):
        """
        """
        r = DnaSequenceAnnotation.objects.filter( subComponent=self.id)
        s = r.count()

        return s     

    def number_related_ParentChildAnnotations( self ):
        """
        """            
        r = DnaSequenceAnnotation.objects.filter( subComponent=self.id)
        s = r.count()
        g = DnaSequenceAnnotation.objects.filter( componentAnnotated=self.id)
        total = s+g.count()
    
        return total  

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

  

class ChemicalComponent(Component):
    """
    Description of a chemical.
    """
    #: To define
    componentType = models.ManyToManyField(ChemicalComponentType, 
                                           blank=True, null=True, 
                                           verbose_name='Part type', 
                                           help_text='Classification of this part.')   
    variantOf = models.ManyToManyField( 'self', symmetrical=False, 
                                            blank=True, null=True, 
                                            verbose_name='Variant of', 
                                            help_text='Specify part(s) this part is derived from, if any.' )
        
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'chemicalcomponent/%i/' % self.id

    class Meta:
        app_label = 'labrack'        
        verbose_name = 'Chemical'


class Chassis(Component):
    """
    Description of a host system. Usually this will be a cell type or bacterial
    strain.
    """
    componentType = models.ManyToManyField(ChassisComponentType, 
                                           blank=True, null=True, 
                                           verbose_name='Type', 
                                           help_text='Classification of this part.')
    variantOf = models.ManyToManyField( 'self', symmetrical=False, 
                                            blank=True, null=True, 
                                            verbose_name='Variant of', 
                                            help_text='Specify part(s) this part is derived from, if any.' )
        
    
    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'chassis/%i/' % self.id
    
    def get_absolute_url(self):
        """requires site URL to be set in admin/sites web interface"""
        return '/chassis/%s/' % self.displayId

    def relatedChassisSamples( self ):
        """
        """
        import labrack.models as M
        r = M.ChassisSample.objects.filter(chassis=self.id)
        return r
    def relatedDnaSamples(self):
        """
        """       
        import labrack.models as M
        r = M.DnaSample.objects.filter(inChassis=self.id)

        return r       
    
    class Meta:
        app_label = 'labrack'           
        verbose_name = 'Cell'



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

    class Meta:
        app_label = 'labrack'    

        
class PeptideComponent(Component):
    """
    Description of a peptide 'part'.
    """
    sequence = models.TextField( help_text='amino acid sequence', 
                                 blank=True, null=True )

    componentType = models.ManyToManyField(PeptideComponentType, 
                                           blank=True, null=True, 
                                           verbose_name='Part type', 
                                           help_text='Classification of this part.')   
    variantOf = models.ManyToManyField( 'self', symmetrical=False, 
                                            blank=True, null=True, 
                                            verbose_name='Variant of', 
                                            help_text='Specify part(s) this part is derived from, if any.' )

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

    class Meta:
        app_label = 'labrack'      
        verbose_name = 'Peptide'


class ProteinComponent(Component):
    """
    Description of a protein 'part'.
    """
    #: optional sequence
    sequence = models.TextField( help_text='amino acid sequence', 
                                blank=True, null=True )
    
    componentType = models.ManyToManyField(ProteinComponentType, 
                                            blank=True, null=True, 
                                            verbose_name='Part type', 
                                            help_text='Classification of this part.')   
    
    GenBankfile = models.FileField(upload_to='documents/GenBank/%Y/%m/%d',blank=True,null=True)
    
    variantOf = models.ManyToManyField( 'self', symmetrical=False, 
                                            blank=True, null=True, 
                                            verbose_name='Variant of', 
                                            help_text='Specify part(s) this part is derived from, if any.' )
  
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
    
    def save(self, *args, **kwargs):
        #Saving the sequence
        self.sequence = self.related_seq()
        super(ProteinComponent, self).save(*args, **kwargs) # Call the "real" save() method.
        if self.GenBankfile:
            self.save_annotation()
    
    def saveWithoutAnnotations(self, *args, **kwargs):
        super(ProteinComponent, self).save(*args, **kwargs) # Call the "real" save() method.    
    
    class Meta :
        app_label = 'labrack'                   
        verbose_name = 'Protein part'
        
class DnaComponent(Component):
    """
    Description of a stretch of DNA or RNA.
    """
   
    sequence = models.TextField( help_text='nucleotide sequence', 
                                             blank=True, null=True )      

    translatesTo = models.ForeignKey( 'ProteinComponent', blank=True, null=True, 
                                      related_name='encodedBy', 
                                      help_text='Protein part this sequence translates to' )

    optimizedFor = models.ForeignKey( 'Chassis', blank=True, null=True )

    componentType = models.ManyToManyField(DnaComponentType, 
                                           blank=True, null=True, 
                                           verbose_name='Type', 
                                           help_text='Classification of this part.')

    circular = models.BooleanField( 'Circular', default=True, 
                                    help_text='is the DNA Circular or not.')


    GenBankfile = models.FileField(upload_to='documents/GenBank/%Y/%m/%d',blank=True,null=True)
    
    
    variantOf = models.ManyToManyField( 'self', symmetrical=False, 
                                            blank=True, null=True, 
                                            verbose_name='Variant of', 
                                            help_text='Specify part(s) this part is derived from, if any.' )
    
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

    def show_resistance(self):
        """filter selection marker for this dnaComponent
        """
        ret = ""
        r = DnaSequenceAnnotation.objects.filter( subComponent=self.id).order_by('bioStart')
                      
        for s in r:
            com = s.componentAnnotated
            comDna= DnaComponent.objects.get(id = com.id)
            for t in comDna.componentType.all():
                if t.name == 'SelectionMarker':
                    if ret=="":
                        ret += s.componentAnnotated.name
                    else:
                        ret += ","+s.componentAnnotated.name
     

        return ret   
    show_resistance.short_description = 'Resistance'
    
    def getVector(self):
        """
        """
        ret = ""
        r = DnaSequenceAnnotation.objects.filter( subComponent=self.id).order_by('bioStart')
        for s in r:
            com = s.componentAnnotated
            comDna= DnaComponent.objects.get(id = com.id)
            for t in comDna.componentType.all():
                if t.name == 'Vector':
                    return s.componentAnnotated
 

        return ret   
    getVector.short_description = 'Vector'    
    
    def related_dnaSamples(self):
        """
        """       
        import labrack.models as M

        r = M.DnaSample.objects.filter(dnaConstruct=self.id)

        return r
    show_optimizedFor.short_description = 'optimized for'
    show_optimizedFor.admin_order_field = 'optimizedFor'

        
    def related_seq( self ):
        """
        """
        gb_features = ''

        gb_file = settings.MEDIA_ROOT+"/"+os.path.normpath(self.GenBankfile.name)
        for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
            gb_features += gb_record.seq.tostring()
        return gb_features
 

    def saveSequenceWithoutAnnotations(self, *args, **kwargs):
        #Saving the sequence 
        super(DnaComponent, self).save(*args, **kwargs) # Call the "real" save() method.
        if self.GenBankfile:
            self.sequence = self.related_seq()
        super(DnaComponent, self).save(*args, **kwargs) # Call the "real" save() method.


    def saveWithoutAnnotations(self, *args, **kwargs):
        
        super(DnaComponent, self).save(*args, **kwargs) # Call the "real" save() method.

    def save_annotation( self ):
        if (self.GenBankfile):
            gb_file = settings.MEDIA_ROOT+"/"+os.path.normpath(self.GenBankfile.name)
            gb_features = ""
            dispId = 1
            isParsingDone = False

            for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
                if (not isParsingDone):
                    for ind in xrange(len(gb_record.features)) :
                        isParsingDone = True
                        nameType = repr(gb_record.features[ind].type).replace("'", "")
                        strandValue = repr(gb_record.features[ind].strand)
                        startPos = repr(gb_record.features[ind].location._start.position+1)
                        endPos = repr(gb_record.features[ind].location._end.position)
                        label = repr(gb_record.features[ind].qualifiers.get('label')).replace("['","").replace("']","").replace("\\"," ")
                        if (label == 'None' ):
                            label = repr(gb_record.features[ind].qualifiers.get('gene')).replace("['","").replace("']","").replace("\\"," ")

                        ###check if the annotation refering to this dna already exists using start and end point and dnaID (this is to solve a bug but should be removed, the bug is that its saving the annotated dna twice))
                        if (not DnaSequenceAnnotation.objects.filter(bioStart = startPos, bioEnd = endPos, strand = strandValue, subComponent = self)):
                            fullSequence = gb_record.seq.tostring()
                            if (startPos==endPos):
                                partOfSequence = fullSequence[int(startPos):int(endPos)+1].replace(" ","").upper()
                            else:
                                partOfSequence = fullSequence[int(startPos):int(endPos)].replace(" ","").upper()

                            if (strandValue == '-1'):
                                partOfSequence = Seq(partOfSequence).reverse_complement().tostring()

                            # the reason to save it twice is to get a unique ID to be able to put it in DisplayId
                            # retrieve the part type, if not existing create it
                            if (not DnaComponentType.objects.filter(name='GenebankType')):
                                subCtGenebankType = DnaComponentType(name = 'GenebankType')
                                subCtGenebankType.save() 
                            subCtGenebankType = DnaComponentType.objects.get(name='GenebankType')
                            if (not DnaComponentType.objects.filter(name=nameType)):
                                ct = DnaComponentType(name = nameType)
                                ct.save()
                            ct = DnaComponentType.objects.filter(name=nameType)
                            ct.subTypeOf = subCtGenebankType

                            if (not self.__class__.objects.filter(sequence=partOfSequence)):
                                dna2db = self.__class__(displayId='########', sequence = partOfSequence, name = label, GenBankfile = None)
                                dna2db.saveWithoutAnnotations()
                                dna2db.componentType = ct
                                dna2db.displayId="gb%06i"%dna2db.id
                                dna2db.saveWithoutAnnotations()
                            dna2db = self.__class__.objects.get(sequence=partOfSequence)
                            # save the annotation in the database
                            an2db = DnaSequenceAnnotation(uri ='', bioStart = startPos, bioEnd = endPos, strand = strandValue, subComponent = self, componentAnnotated = dna2db)
                            an2db.save()             
    class Meta:
        app_label = 'labrack'
        verbose_name = 'DNA part'



