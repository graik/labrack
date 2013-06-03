from django.db import models
from usermixin import UserMixinModel
from datetime import datetime

# Create your models here.

class Location(UserMixinModel):
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
        
          
        
        



class Rack(UserMixinModel):
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

        return r    

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'rack/%i/' % self.id          
    class Meta:
        app_label = 'labrack'   


class Container( UserMixinModel ):
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

    def __unicode__( self ):
        return "%s  ,%s, %s" % (self.displayId,self.rack,self.rack.current_location)

    def get_relative_url(self):
        """
        Define standard relative URL for object access in templates
        """
        return 'container/%i/' % self.id    

    class Meta:
        app_label = 'labrack'   
        ordering = ('displayId',)





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
    
    class Meta:
            app_label = 'labrack'      


    