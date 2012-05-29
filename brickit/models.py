from django.db import models

from django.contrib.auth.models import User

# Create your models here.
APP_URL = '/registry'

# Create your models here.
class StorageContainer( models.Model ):
    """
    A container holding several physical samples of DNA or other stuff.
    """

    STORAGE_CONTAINER_TYPES= (
        ('96-well-plate', '96 well plate'),
        ('384-well-plate','384 well plate'),
        ('box', 'freezer box'),
        ('other', 'other' ) )

    #: actual label stuck to this container
    label = models.CharField('label', max_length=20, unique=True,
                             help_text='example: D_001 or C_012 (D...DNA, C...Cells)')

    #: [optional] annotation
    description = models.CharField(max_length=200, blank=True )

    #: [optional]
    comments = models.TextField(blank=True)

    #: what type of container is it
    container_type = models.CharField('type of container', max_length=30,
                                      choices=STORAGE_CONTAINER_TYPES )

    location = models.CharField('location', max_length=200)

    #: creators or users
    users = models.ManyToManyField(User, null=True, blank=True, db_index=True)


    def __unicode__( self ):
        return "%s" % self.label


    def get_absolute_url(self):
        """Define standard URL for object.get_absolute_url access in templates """
        return APP_URL+'/storagecontainer/%i/' % self.id
 
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
        ordering = ('label',)

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
    label = models.CharField('label', max_length=20,
                             help_text='example: rg23/12/07-1a (rg...user initials)' )

    #: link to a single container
    container = models.ForeignKey( StorageContainer, related_name='samples' )

    well_type = models.CharField('type of well or tube', max_length=30,
                                 choices=STORAGE_WELL_TYPES,
                                 default='tube')

    storage_type = models.CharField('stored as', max_length=100,
                                    choices=STORAGE_TYPES, default='dna' )

    #: creators or users
    users = models.ManyToManyField(User, db_index=True)

    #: [optional]
    comments = models.TextField(blank=True)

    #: [optional] (for now)
    well = models.CharField(max_length=5, blank=True )

    created = models.DateField('created at', auto_now_add=True)

    #: link to the physical DNA contained in this sample
    dna = models.ForeignKey( 'DNA',
                             verbose_name='physical DNA',
                             related_name='in_samples',
                             blank=False)

    #: [optional] strain this DNA is stored in, if any
    cell = models.ForeignKey( 'Cell',
                              verbose_name='in cell',
                              related_name='sample',
                              blank=True,
                              null=True)

    #: concentration in ng / ul (=mg/l)
    concentration = models.DecimalField( max_digits=6, decimal_places=2,
                                         default=50 )
    #: unit of given concentration
    concentration_unit = models.CharField( 'conc. unit', max_length=10,
                                           choices=CONCENTRATION_UNITS,
                                           default='mg/l')

    def delete( self ):
        """
        Delete this sample but, if possible, also the connected Physical DNA.
        Note: this function is called but the delete somehow doesn't work.
        """
        dna = self.dna
        super( Sample, self).delete()
        if dna:
            if not dna.in_samples: ## no other samples refer to this DNA
                dna.delete()

    def __unicode__( self ):
        return self.container.label + ' / ' + self.label

    def get_absolute_url(self):
        """Define standard URL for object.get_absolute_url access in templates """
        return APP_URL+'/sample/%i/' % self.id

    def get_sequencing_evaluation( self ):
        """
        @return: str, best registered sequencing analysis 
                ('confirmed', 'inconsistent', 'ambiguous' or 'not analyzed')
        """
        results = [ seq.evaluation for seq in self.sequencing.all() ]

        for verdict in Sequencing.EVALUATION:
            if verdict[0] in results:
                return verdict[1]

        return 'none'
    get_sequencing_evaluation.short_description = 'Sequencing'

    def related_samples( self ):
        """
        @return: QuerySet of samples with similar DNA content
        """
        if self.dna.biobrick:
            r = Sample.objects.filter( dna__biobrick=self.dna.biobrick )
        else:
            r = Sample.objects.filter( dna__biobrick=self.dna.biobrick,
                                       dna__vector=self.dna.vector )
            
        return r.exclude( pk=self.pk )
    
    def get_sequence( self, recenter=0 ):
        return self.dna.get_sequence( recenter=recenter )

    def get_pretty_sequence( self, recenter=0 ):
        return self.dna.get_pretty_sequence( recenter=recenter )

    class Meta:
        unique_together = ('label', 'container')
        ordering = ('container', 'well')

