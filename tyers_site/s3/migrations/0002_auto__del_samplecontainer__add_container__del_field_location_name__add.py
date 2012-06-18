# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting model 'SampleContainer'
        db.delete_table('s3_samplecontainer')

        # Removing M2M table for field users on 'SampleContainer'
        db.delete_table('s3_samplecontainer_users')

        # Adding model 'Container'
        db.create_table('s3_container', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('containerType', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('location', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['s3.Location'])),
            ('owner', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='owner', null=True, to=orm['auth.User'])),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('s3', ['Container'])

        # Adding M2M table for field users_read on 'Container'
        db.create_table('s3_container_users_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['s3.container'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('s3_container_users_read', ['container_id', 'user_id'])

        # Adding M2M table for field users_write on 'Container'
        db.create_table('s3_container_users_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['s3.container'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('s3_container_users_write', ['container_id', 'user_id'])

        # Adding M2M table for field group_read on 'Container'
        db.create_table('s3_container_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['s3.container'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('s3_container_group_read', ['container_id', 'group_id'])

        # Adding M2M table for field group_write on 'Container'
        db.create_table('s3_container_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['s3.container'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('s3_container_group_write', ['container_id', 'group_id'])

        # Deleting field 'Location.name'
        db.delete_column('s3_location', 'name')

        # Adding field 'Location.displayId'
        db.add_column('s3_location', 'displayId',
                      self.gf('django.db.models.fields.CharField')(default=1, unique=True, max_length=20),
                      keep_default=False)

        # Adding field 'Location.shortDescription'
        db.add_column('s3_location', 'shortDescription',
                      self.gf('django.db.models.fields.CharField')(default='', max_length=200, blank=True),
                      keep_default=False)


        # Changing field 'Location.room'
        db.alter_column('s3_location', 'room', self.gf('django.db.models.fields.CharField')(max_length=20))

        # Changing field 'Sample.container'
        db.alter_column('s3_sample', 'container_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['s3.Container']))

    def backwards(self, orm):
        # Adding model 'SampleContainer'
        db.create_table('s3_samplecontainer', (
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(max_length=20, unique=True)),
            ('containerType', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('s3', ['SampleContainer'])

        # Adding M2M table for field users on 'SampleContainer'
        db.create_table('s3_samplecontainer_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('samplecontainer', models.ForeignKey(orm['s3.samplecontainer'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('s3_samplecontainer_users', ['samplecontainer_id', 'user_id'])

        # Deleting model 'Container'
        db.delete_table('s3_container')

        # Removing M2M table for field users_read on 'Container'
        db.delete_table('s3_container_users_read')

        # Removing M2M table for field users_write on 'Container'
        db.delete_table('s3_container_users_write')

        # Removing M2M table for field group_read on 'Container'
        db.delete_table('s3_container_group_read')

        # Removing M2M table for field group_write on 'Container'
        db.delete_table('s3_container_group_write')


        # User chose to not deal with backwards NULL issues for 'Location.name'
        raise RuntimeError("Cannot reverse this migration. 'Location.name' and its values cannot be restored.")
        # Deleting field 'Location.displayId'
        db.delete_column('s3_location', 'displayId')

        # Deleting field 'Location.shortDescription'
        db.delete_column('s3_location', 'shortDescription')


        # Changing field 'Location.room'
        db.alter_column('s3_location', 'room', self.gf('django.db.models.fields.CharField')(max_length=25))

        # Changing field 'Sample.container'
        db.alter_column('s3_sample', 'container_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['s3.SampleContainer']))

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
        's3.chassis': {
            'Meta': {'object_name': 'Chassis'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.Chassis']", 'null': 'True', 'blank': 'True'})
        },
        's3.collection': {
            'Meta': {'object_name': 'Collection'},
            'chassis': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.Chassis']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'dnaComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.DnaComponent']", 'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'proteinComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.ProteinComponent']", 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        's3.componenttype': {
            'Meta': {'object_name': 'ComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['s3.ComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200'})
        },
        's3.container': {
            'Meta': {'ordering': "('displayId',)", 'object_name': 'Container'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'groups_read'", 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'groups_write'", 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'location': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['s3.Location']"}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'owner'", 'null': 'True', 'to': "orm['auth.User']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'users_read': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'users_read'", 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'users_write': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'users_write'", 'symmetrical': 'False', 'to': "orm['auth.User']"})
        },
        's3.dnacomponent': {
            'Meta': {'object_name': 'DnaComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['s3.Chassis']", 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['s3.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.DnaComponent']", 'null': 'True', 'blank': 'True'})
        },
        's3.location': {
            'Meta': {'object_name': 'Location'},
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        's3.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.ProteinComponent']", 'null': 'True', 'blank': 'True'})
        },
        's3.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chassis_samples'", 'null': 'True', 'to': "orm['s3.Chassis']"}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '50', 'max_digits': '6', 'decimal_places': '2'}),
            'concentrationUnit': ('django.db.models.fields.CharField', [], {'default': "'mg/l'", 'max_length': '10'}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['s3.Container']"}),
            'created': ('django.db.models.fields.DateField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'dna': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'dna_samples'", 'null': 'True', 'to': "orm['s3.DnaComponent']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'protein_samples'", 'null': 'True', 'to': "orm['s3.ProteinComponent']"}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'vector': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'vector_samples'", 'null': 'True', 'to': "orm['s3.VectorDnaComponent']"}),
            'vesselType': ('django.db.models.fields.CharField', [], {'default': 'None', 'max_length': '30', 'null': 'True', 'blank': 'True'})
        },
        's3.selectivemarker': {
            'Meta': {'ordering': "('name',)", 'object_name': 'SelectiveMarker'},
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        },
        's3.sequenceannotation': {
            'Meta': {'object_name': 'SequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'precedes': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1', 'null': 'True', 'blank': 'True'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        's3.vectordnacomponent': {
            'Meta': {'object_name': 'VectorDnaComponent', '_ormbases': ['s3.DnaComponent']},
            'dnacomponent_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['s3.DnaComponent']", 'unique': 'True', 'primary_key': 'True'}),
            'marker': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['s3.SelectiveMarker']", 'null': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['s3']