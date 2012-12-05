from django.conf.urls.defaults import patterns, url
from labrack.views import dnalist,list
 

urlpatterns = patterns('labrack.views',
    url(r'dnalist/$', dnalist, name='dnalist'),
    url(r'list/$', list, name='list') 
)
