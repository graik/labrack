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


import os.path as osp

def absfile( filename, resolveLinks=1 ):
    """
    Get absolute file path::
      - expand ~ to user home, change
      - expand ../../ to absolute path
      - resolve links
      - add working directory to unbound files ('ab.txt'->'/home/raik/ab.txt')

    @param filename: name of file
    @type  filename: str
    @param resolveLinks: eliminate any symbolic links (default: 1)
    @type  resolveLinks: 1|0
    
    @return: absolute path or filename
    @rtype: string

    @raise ToolsError: if a ~user part does not translate to an existing path
    """
    if not filename:
        return filename
    r = osp.abspath( osp.expanduser( filename ) )

    if '~' in r:
        raise ToolsError, 'Could not expand user home in %s' % filename

    if resolveLinks:
        r = osp.realpath( r )
    r = osp.normpath(r)
    return r


def approot():
    """
    Root of labrack project.
    
    @return: absolute path of the root of current project::
             i.e. '/home/Bis/raik/biskit'
    @rtype: string
    """
    ## import this module
    #from labrack import tools
    ## get location of this module
    f = absfile(__file__)
    ## extract path and assume it is 'project_root/Biskit'
    f = osp.split( f )[0] + '/../'
    return absfile( f )
