from django.conf.urls.defaults import patterns, include, url

from django.contrib import admin

from labrack.jsViews import getGBDataFromFile,fileUpload,\
     getDnaInfo,searchDnaParts,getAnnotToBeDeleted

#import labrack.componentViews as CViews
import settings

# deprecated in django 1.5: Use class-based views instead
# http://stackoverflow.com/questions/11428427/no-module-named-simple-error-in-django
# from django.views.generic.simple import redirect_to
admin.autodiscover()

urlpatterns = patterns('',
    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    
    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
    #url(r'^labrack/chassis[/]?$', CViews.chassis_list),
    #url(r'^labrack/chassis/(?P<id>[0-9A-Za-z]+)', CViews.ChassisDetail.as_view()),
    #url(r'^labrack/dna/(?P<id>[0-9A-Za-z]+)', CViews.DnaComponentDetail.as_view()),
    url( r'fileUpload/$', fileUpload, name="fileUpload" ),
    url(r'^getDnaInfo/(?P<text_value>.*)/$',getDnaInfo,name='getDnaInfo'),
    url(r'^getGBDataFromFile/(?P<filePath>.*)/$', getGBDataFromFile,name='getGBDataFromFile'),  
    url(r'^searchDnaParts/(?P<sequence_text>.*)/(?P<displayIdDnaComponent>.*)/$', searchDnaParts,name='searchDnaParts'),  
    url(r'^checkAnnotToDelete/(?P<jsonAmmpt>.*)/(?P<displayIdDnaComponent>.*)/$', getAnnotToBeDeleted,name='getAnnotToBeDeleted'),      
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT})
)

