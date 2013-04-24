## Copyright 2012-2013 Raik Gruenberg

## This file is part of the labrack project (http://labrack.sf.net)
## labrack is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.

## labrack is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public
## License along with labrack. If not, see <http://www.gnu.org/licenses/>.

"""Views for display (rather than editing) of Components

A perhaps better short term solution is this:
http://stackoverflow.com/questions/6680631/django-admin-separate-read-only-view-and-change-view
This would retain everything within Admin and create a view form separate
from the change/edit form.
"""


from django.views.generic import ListView, DetailView

import labrack.models as M
import labrack.admin as A
from labrack.displayAdmin import labracksite



class ComponentDetailView( DetailView ):
    """
    Class-based view for Component pages.
    
    The default DetailView expects that the primary key (pk) is part of the
    request URL and that objects are fetched based on this primary key. We
    instead want to reach components by their displayID, e.g.:
    labrack/dnacomponent/sb0010 instead of
    admin/labrack/dnacomponent/2
    That's why ComponentDetailView overrides the method for fetching objects.
    """

    def get(self, request, *args, **kwargs):
        """Override BaseDetailView.get to use displayID instead of pk"""
        
        self.object = self.model.objects.get(displayId=kwargs['id'])
        context = self.get_context_data(object=self.object)
        return self.render_to_response(context)
    
    
class ChassisDetail( ComponentDetailView ):
    model = M.Chassis
    
class DnaComponentDetail( ComponentDetailView ):
    model = M.DnaComponent
    
    

def chassis_list( request, *args, **kwargs):
    """
    View capturing labrack/chassis and returning a ChangeList View from
    a customized adminSite instance (labracksite).
    
    Kind of works, but all links still point to /admin/labrack/chassis rather
    than /labrack/chassis.
    """
    a = A.ChassisAdmin( M.Chassis, labracksite )
    return a.changelist_view( request, extra_context=None )