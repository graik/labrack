# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding M2M table for field owners on 'Sample'
        db.create_table('labrack_sample_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sample', models.ForeignKey(orm['labrack.sample'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_sample_owners', ['sample_id', 'user_id'])

        # Adding M2M table for field group_read on 'Sample'
        db.create_table('labrack_sample_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sample', models.ForeignKey(orm['labrack.sample'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_sample_group_read', ['sample_id', 'group_id'])

        # Adding M2M table for field group_write on 'Sample'
        db.create_table('labrack_sample_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sample', models.ForeignKey(orm['labrack.sample'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_sample_group_write', ['sample_id', 'group_id'])


    def backwards(self, orm):
        # Removing M2M table for field owners on 'Sample'
        db.delete_table('labrack_sample_owners')

        # Removing M2M table for field group_read on 'Sample'
        db.delete_table('labrack_sample_group_read')

        # Removing M2M table for field group_write on 'Sample'
        db.delete_table('labrack_sample_group_write')


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
        'labrack.chassis': {
            'Meta': {'object_name': 'Chassis', '_ormbases': ['labrack.Component']},
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.chemicalcomponent': {
            'Meta': {'object_name': 'ChemicalComponent', '_ormbases': ['labrack.Component']},
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.collection': {
            'Meta': {'object_name': 'Collection'},
            'components': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Component']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.component': {
            'Meta': {'object_name': 'Component'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'component_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'component_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'component_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'component_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Component']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.componenttype': {
            'Meta': {'object_name': 'ComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.container': {
            'Meta': {'ordering': "('displayId',)", 'object_name': 'Container'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'container_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'location': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'containers'", 'to': "orm['labrack.Location']"}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"})
        },
        'labrack.dnacomponent': {
            'Meta': {'object_name': 'DnaComponent', '_ormbases': ['labrack.Component']},
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"})
        },
        'labrack.location': {
            'Meta': {'object_name': 'Location'},
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        'labrack.peptidecomponent': {
            'Meta': {'object_name': 'PeptideComponent', '_ormbases': ['labrack.Component']},
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent', '_ormbases': ['labrack.Component']},
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'})
        },
        'labrack.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'aliquotNr': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'attachment': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['labrack.Container']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sample_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'preparation_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2012, 7, 20, 0, 0)'}),
            'sampleCollection': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SampleCollection']", 'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'ok'", 'max_length': '30'})
        },
        'labrack.samplecollection': {
            'Meta': {'object_name': 'SampleCollection'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplecollection_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'samplecollection_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'samplecollection_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '25'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'samplecollection_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"})
        },
        'labrack.samplecontent': {
            'Meta': {'object_name': 'SampleContent'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amountUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'amountUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentrationUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'concentrationUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sampleContent'", 'to': "orm['labrack.Sample']"})
        },
        'labrack.samplelink': {
            'Meta': {'object_name': 'SampleLink'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'link': ('django.db.models.fields.CharField', [], {'max_length': '1000'}),
            'linkType': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samplelink'", 'to': "orm['labrack.Sample']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'})
        },
        'labrack.sampleprovenance': {
            'Meta': {'object_name': 'SampleProvenance'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'provenanceType': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sampleProvenance'", 'to': "orm['labrack.Sample']"}),
            'sample_source': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sampleSource'", 'null': 'True', 'to': "orm['labrack.Sample']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'})
        },
        'labrack.selectivemarker': {
            'Meta': {'ordering': "('name',)", 'object_name': 'SelectiveMarker'},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'})
        },
        'labrack.sequenceannotation': {
            'Meta': {'object_name': 'SequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'precedes': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1', 'null': 'True', 'blank': 'True'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.unit': {
            'Meta': {'object_name': 'Unit'},
            'conversion': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'unitType': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        }
    }

    complete_apps = ['labrack']