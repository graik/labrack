from django.conf.urls.defaults import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()
import settings

urlpatterns = patterns('',
    # Examples:

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    #url(r'^labrack/', include('labrack.urls')),
    #url(r'^admin/labrack/dnacomponent/(?P<Id>\d+)/$', index, name='index'),    
    url(r'^', include('labrack.urls')),    
    url(r'^admin/', include(admin.site.urls)),
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),
    #url(r'^sampleFiles/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.FOLDER_FILES_PATH + '/sampleFiles'}),
    #url(r'sequenceannotation/add/$', SequenceAnnotationCreate.as_view(), name='sequenceannotation_add'),
    #url(r'sequenceannotation/(?P<pk>\d+)/$', SequenceAnnotationUpdate.as_view(), name='sequenceannotation_update'),
    #url(r'sequenceannotation/(?P<pk>\d+)/delete/$', SequenceAnnotationDelete.as_view(), name='sequenceannotation_delete'),
    #url(r'^documents/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_URL+ '/documents'}),
    #url(r'^admin/labrack/$', index, name='index')
)
"""
urlpatterns = patterns('',
    url(r'^$', views.index, name='index')
)
"""