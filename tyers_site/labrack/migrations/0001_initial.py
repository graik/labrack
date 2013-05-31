# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Location'
        db.create_table('labrack_location', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='location_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('temperature', self.gf('django.db.models.fields.FloatField')()),
            ('room', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 5, 27, 0, 0))),
        ))
        db.send_create_signal('labrack', ['Location'])

        # Adding M2M table for field owners on 'Location'
        db.create_table('labrack_location_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('location', models.ForeignKey(orm['labrack.location'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_location_owners', ['location_id', 'user_id'])

        # Adding model 'Rack'
        db.create_table('labrack_rack', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='rack_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('current_location', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Location'], null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 5, 27, 0, 0))),
        ))
        db.send_create_signal('labrack', ['Rack'])

        # Adding M2M table for field owners on 'Rack'
        db.create_table('labrack_rack_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('rack', models.ForeignKey(orm['labrack.rack'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_rack_owners', ['rack_id', 'user_id'])

        # Adding model 'Container'
        db.create_table('labrack_container', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='container_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('containerType', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('rack', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Rack'])),
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 5, 27, 0, 0))),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['Container'])

        # Adding M2M table for field owners on 'Container'
        db.create_table('labrack_container_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['labrack.container'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_container_owners', ['container_id', 'user_id'])


    def backwards(self, orm):
        # Deleting model 'Location'
        db.delete_table('labrack_location')

        # Removing M2M table for field owners on 'Location'
        db.delete_table('labrack_location_owners')

        # Deleting model 'Rack'
        db.delete_table('labrack_rack')

        # Removing M2M table for field owners on 'Rack'
        db.delete_table('labrack_rack_owners')

        # Deleting model 'Container'
        db.delete_table('labrack_container')

        # Removing M2M table for field owners on 'Container'
        db.delete_table('labrack_container_owners')


    models = {
        'auth.group': {
            'Meta': {'object_name': 'Group'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        'auth.permission': {
            'Meta': {'ordering': "('content_type__app_label', 'content_type__model', 'codename')", 'unique_together': "(('content_type', 'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        'labrack.container': {
            'Meta': {'ordering': "('displayId',)", 'object_name': 'Container'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'container_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'rack': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Rack']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 5, 27, 0, 0)'})
        },
        'labrack.location': {
            'Meta': {'object_name': 'Location'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'location_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'location_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 5, 27, 0, 0)'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        'labrack.rack': {
            'Meta': {'object_name': 'Rack'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'rack_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'current_location': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Location']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'rack_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 5, 27, 0, 0)'})
        }
    }

    complete_apps = ['labrack']