from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
from labrack.views import getTypeDnaInfo,reviewdna,getParentTypeDnaInfo
from django.contrib.auth.decorators import login_required
from django.contrib import databrowse
from labrack.search.search_view import basic_search
import settings






admin.autodiscover()

urlpatterns = patterns('',
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^getTypeDnaInfo/(?P<maintype>.*)/$',getTypeDnaInfo,name='getTypeDnaInfo'),
    url(r'^getParentTypeDnaInfo/(?P<subtype>.*)/$',getParentTypeDnaInfo,name='getParentTypeDnaInfo'),
    #url(r'^search/', include('haystack.urls')),
    url(r'^search/', basic_search,name='basic_search'),
    url(r'^reviewdna/(?P<displayId>.*)/$',reviewdna,name='reviewdna'),
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT})    
    )