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

from django.contrib.admin import SimpleListFilter
from django.db.models import Q
from django.utils.translation import ugettext_lazy as _



class ContainerListFilter(SimpleListFilter):
    # Human-readable title which will be displayed in the
    # right admin sidebar just above the filter options.
    title = _('container')

    # Parameter for the filter that will be used in the URL query.
    parameter_name = 'container__id__exact'

    def lookups(self, request, model_admin):
        """
        Returns a list of tuples. The first element in each
        tuple is the coded value for the option that will
        appear in the URL query. The second element is the
        human-readable name for the option that will appear
        in the right sidebar.
        """
        
        qs = model_admin.queryset(request)
        
        containerList = []
        containerDict = {} 
        
        for sample in qs:
            if containerDict.has_key(sample.container.pk) != True:
                containerList.append(( str(sample.container.pk), _(str(sample.container)) ))
            containerDict[sample.container.pk] = True

        return containerList


    def queryset(self, request, qs):
        
        if request.user.is_superuser:
            pass
        else:
            pass
        
        if(self.value() == None):
            return qs
        else:
            return qs.filter(container__id__exact=self.value())





class SampleCollectionListFilter(SimpleListFilter):
    # Human-readable title which will be displayed in the
    # right admin sidebar just above the filter options.
    title = _('sampleCollection')

    # Parameter for the filter that will be used in the URL query.
    parameter_name = 'sampleCollection__id__exact'

    def lookups(self, request, model_admin):
        """
        Returns a list of tuples. The first element in each
        tuple is the coded value for the option that will
        appear in the URL query. The second element is the
        human-readable name for the option that will appear
        in the right sidebar.
        """
        
        qs = model_admin.queryset(request)
        
        sampleCollectionList = []
        sampleCollectionDict = {} 
        
        for sample in qs:
            for sampleCollection in sample.sampleCollection.all():
                if sampleCollectionDict.has_key(sampleCollection.pk) != True:
                    sampleCollectionList.append(( str(sampleCollection.pk), _(str(sampleCollection)) ))
                    sampleCollectionDict[sampleCollection.pk] = True

        return sampleCollectionList


    def queryset(self, request, qs):
        
        if request.user.is_superuser:
            pass
        else:
            qs = qs.filter(Q(owners=request.user) |
                           Q(group_read=request.user.groups.all()) |
                           Q(group_write=request.user.groups.all()) 
                           ).distinct()
        
        if(self.value() == None):
            return qs
        else:
            return qs.filter(sampleCollection__id__exact=self.value())