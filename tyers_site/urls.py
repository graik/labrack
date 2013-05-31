from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
from labrack.views import getTypeDnaInfo,reviewdna,getParentTypeDnaInfo
from django.contrib.auth.decorators import login_required

from django.contrib import databrowse


admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'tyers_site.views.home', name='home'),
    # url(r'^tyers_site/', include('tyers_site.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^getTypeDnaInfo/(?P<maintype>.*)/$',getTypeDnaInfo,name='getTypeDnaInfo'),
    url(r'^getParentTypeDnaInfo/(?P<subtype>.*)/$',getParentTypeDnaInfo,name='getParentTypeDnaInfo'),
    url(r'^search/', include('haystack.urls')),
    url(r'^reviewdna/(?P<displayId>.*)/$',reviewdna,name='reviewdna')
    )