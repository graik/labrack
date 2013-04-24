from django.conf.urls.defaults import patterns, include, url

from django.contrib import admin

from labrack.views import getAnnotToBeDeleted,getGBDataFromFile,file_upload,\
     get_dna_info,dnalist,celllist,cellonlylist,success,\
     search_dna_parts

import labrack.componentViews as CViews
import settings

# deprecated in django 1.5: Use class-based views instead
# http://stackoverflow.com/questions/11428427/no-module-named-simple-error-in-django
# from django.views.generic.simple import redirect_to


urlpatterns = patterns('',
    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    
    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
    url(r'^labrack/chassis[/]?$', CViews.chassis_list),
    url(r'^labrack/chassis/(?P<id>[0-9A-Za-z]+)', CViews.ChassisDetail.as_view()),
    url(r'^labrack/dna/(?P<id>[0-9A-Za-z]+)', CViews.DnaComponentDetail.as_view()),
    url( r'file_upload/$', file_upload, name="file_upload" ),
    url(r'^get_dna_info/(?P<text_value>.*)/$',get_dna_info,name='get_dna_info'),
    url(r'^getGBDataFromFile/(?P<filePath>.*)/$', getGBDataFromFile,name='getGBDataFromFile'),  
    url(r'^search_dna_parts/(?P<sequence_text>.*)/(?P<displayIdDnaComponent>.*)/$', search_dna_parts,name='search_dna_parts'),  
    url(r'^checkAnnotToDelete/(?P<jsonAmmpt>.*)/(?P<displayIdDnaComponent>.*)/$', getAnnotToBeDeleted,name='getAnnotToBeDeleted'),      
    url(r'^success/$', success, name='success'),    
    url(r'^cellonlylist/$', cellonlylist, name='cellonlylist'),
    url(r'^dnalist/$', dnalist, name='dnalist'),
    url(r'^celllist/$', celllist, name='celllist'),
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT})
)

