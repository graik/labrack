
from django.db import models


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


    