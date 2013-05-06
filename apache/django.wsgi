import os
import sys
import django.core.handlers.wsgi

sys.path.append('/home/django/py')
os.environ['DJANGO_SETTINGS_MODULE'] = 'tyers_site.settings'

_application = django.core.handlers.wsgi.WSGIHandler()

def application(environ, start_response):
    environ['SCRIPT_NAME'] = '/'
    environ['PATH_INFO'] = environ['SCRIPT_NAME'] + environ['PATH_INFO']
    return _application(environ, start_response)

