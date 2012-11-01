from django.conf.urls.defaults import patterns, url
from tyers_site.labrack.models import Sample

urlpatterns = patterns('labrack.views',
    url(r'list/$', 'list', name='list'),
    #url(r'admin/labrack/sample/$', Sample),
)
