from django.conf.urls.defaults import patterns, url
 

urlpatterns = patterns('labrack.views',
    url(r'list/$', 'list', name='list'),
    #url(r'admin/labrack/sample/$', Sample),
)
