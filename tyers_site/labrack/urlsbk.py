from django.conf.urls.defaults import patterns, url
from labrack.views import dnalist,list,celllist
 

urlpatterns = patterns('labrack.views',
    url(r'^dnalist/$', dnalist, name='dnalist'),
    url(r'^celllist/$', celllist, name='celllist'),
    url(r'list/$', list, name='list') 
)
