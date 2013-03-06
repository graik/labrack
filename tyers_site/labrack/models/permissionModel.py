from django.db import models
from django.contrib.auth.models import User
from django.contrib.auth.models import Group


class PermissionModel( models.Model ):
    #: Permissions

    created_by = models.ForeignKey(User, null=True, blank=True, 
                                   related_name='%(class)s_created_by')

    owners = models.ManyToManyField(User, null=True, blank=True, 
                                    related_name='%(class)s_owners')
    
    creation_date = models.DateTimeField(auto_now_add=True, null=True)
    Ptest_user = models.ForeignKey(User, null=True, blank=True)


    modification_date = models.DateTimeField(auto_now=True, null=True)    
    #group_read = models.ManyToManyField(Group, null=True, blank=True, 
                                        #related_name='%(class)s_groups_read')
    #group_write = models.ManyToManyField(Group, null=True, blank=True, 
                                         #related_name='%(class)s_groups_write')


    #def writePermission(self, currentUser):

        ## Superusers can modify
        #if currentUser.is_superuser:
            #return True

        ## Owners can modify
        #for user in self.owners.all():
            #if user == currentUser:
                #return True

        ## Groups with write can modify
        #for group_obj in self.group_write.all():
            #for group_member in currentUser.groups.all():
                #if group_obj == group_member:
                    #return True        

        #return False

    class Meta:
        app_label = 'labrack'        
        abstract = True
