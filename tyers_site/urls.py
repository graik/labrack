from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
from labrack.views import getTypeDnaInfo,reviewdna,getParentTypeDnaInfo \
     ,getParentTypeCellInfo, getTypeCellInfo, reviewplasmid
from django.contrib.auth.decorators import login_required
import settings






admin.autodiscover()

urlpatterns = patterns('',
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^getTypeDnaInfo/(?P<maintype>.*)/$',getTypeDnaInfo,name='getTypeDnaInfo'),
    url(r'^getTypeCellInfo/(?P<maintype>.*)/$',getTypeCellInfo,name='getTypeCellInfo'),
    url(r'^getParentTypeDnaInfo/(?P<subtype>.*)/$',getParentTypeDnaInfo,name='getParentTypeDnaInfo'),
    url(r'^getParentTypeCellInfo/(?P<subtype>.*)/$',getParentTypeCellInfo,name='getParentTypeCellInfo'),
    url(r'^reviewdna/(?P<displayId>.*)/$',reviewdna,name='reviewdna'),
    url(r'^reviewplasmid/(?P<id>.*)/$',reviewplasmid,name='reviewplasmid'),
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT})    
    )