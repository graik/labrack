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

"""Try to rescue Admin Change Lists for public display"""

from django.contrib.admin.sites import AdminSite

import labrack.models as M
import labrack.admin as A

## secondary Admin site for public views

class LabrackAdmin( AdminSite ):
    """Customized Admin Site version for public interface"""
    pass


labracksite = LabrackAdmin( name='labrack', app_name='labrack')


## register models for public web site
labracksite.register(M.Chassis, A.ChassisAdmin)

