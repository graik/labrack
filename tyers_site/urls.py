from django.conf.urls.defaults import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
from labrack.views import getGBDataFromFile,ajax_upload,upload_page,get_dna_info,dnalist,celllist,cellonlylist,contact,success,plasmidlist,search_dna_parts
admin.autodiscover()
from labrack.models.component import Person 
from django.views.generic.simple import redirect_to
import settings

urlpatterns = patterns('',
    # Examples:

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    #url(r'^labrack/', include('labrack.urls')),
    #url(r'^admin/labrack/dnacomponent/(?P<Id>\d+)/$', index, name='index'), 
    #url(r'^grappelli/', include('grappelli.urls')),
    #url(r'^', include('labrack.urls')),
    url( r'ajax_upload/$', ajax_upload, name="ajax_upload" ),
    url( r'/project/$', upload_page, name="upload_page" ),     
    url(r'^get_dna_info/(?P<text_value>.*)/$',get_dna_info,name='get_dna_info'),
    url(r'^getGBDataFromFile/(?P<filePath>.*)/$',getGBDataFromFile,name='getGBDataFromFile'),  
    url(r'^test/(?P<sequence_text>.*)/$',search_dna_parts,name='search_dna_parts'),  
    url(r'^contact/$', contact, name='contact'),
    url(r'^success/$', success, name='success'),    
    url(r'^cellonlylist/$', cellonlylist, name='cellonlylist'),
    url(r'^plasmidlist/$', plasmidlist, name='plasmidlist'),
    url(r'^dnalist/$', dnalist, name='dnalist'),
    url(r'^celllist/$', celllist, name='celllist'),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT})
    #url(r'^sampleFiles/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.FOLDER_FILES_PATH + '/sampleFiles'}),
    #url(r'sequenceannotation/add/$', SequenceAnnotationCreate.as_view(), name='sequenceannotation_add'),
    #url(r'sequenceannotation/(?P<pk>\d+)/$', SequenceAnnotationUpdate.as_view(), name='sequenceannotation_update'),
    #url(r'sequenceannotation/(?P<pk>\d+)/delete/$', SequenceAnnotationDelete.as_view(), name='sequenceannotation_delete'),
    #url(r'^documents/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_URL+ '/documents'}),
    
)
"""
urlpatterns = patterns('',
    url(r'^$', views.index, name='index')
)
"""