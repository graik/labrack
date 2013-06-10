from django.db import models
from django.contrib.auth.models import User
from django.contrib.auth.models import Group


class UserMixinModel( models.Model ):
    #: Permissions

    createdBy = models.ForeignKey(User, null=True, blank=True, 
                                  related_name='%(class)s_created_by')
    
    class Meta:
        app_label = 'labrack'        
        abstract = True
