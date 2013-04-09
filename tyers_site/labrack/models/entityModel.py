from django.db import models
from django.contrib.auth.models import User
from django.contrib.auth.models import Group


class EntityModel( models.Model ):
    #: Permissions

    created_by = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='%(class)s_created_by')

    owners = models.ManyToManyField(User, null=True, blank=True, 
                                    related_name='%(class)s_owners')
    
        
    

    class Meta:
        app_label = 'labrack'        
        abstract = True
