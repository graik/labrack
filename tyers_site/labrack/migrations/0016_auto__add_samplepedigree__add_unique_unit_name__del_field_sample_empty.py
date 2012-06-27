# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'SamplePedigree'
        db.create_table('labrack_samplepedigree', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample', self.gf('django.db.models.fields.related.ForeignKey')(related_name='sameplepedigree', to=orm['labrack.Sample'])),
            ('amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amount_unit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='samplepedigree_amount_unit', null=True, to=orm['labrack.Unit'])),
            ('concentration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('concentration_unit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='samplepedigree_concentration_unit', null=True, to=orm['labrack.Unit'])),
        ))
        db.send_create_signal('labrack', ['SamplePedigree'])

        # Adding unique constraint on 'Unit', fields ['name']
        db.create_unique('labrack_unit', ['name'])

        # Deleting field 'Sample.empty'
        db.delete_column('labrack_sample', 'empty')

        # Adding field 'Sample.status'
        db.add_column('labrack_sample', 'status',
                      self.gf('django.db.models.fields.CharField')(default='ok', max_length=30),
                      keep_default=False)


    def backwards(self, orm):
        # Removing unique constraint on 'Unit', fields ['name']
        db.delete_unique('labrack_unit', ['name'])

        # Deleting model 'SamplePedigree'
        db.delete_table('labrack_samplepedigree')

        # Adding field 'Sample.empty'
        db.add_column('labrack_sample', 'empty',
                      self.gf('django.db.models.fields.BooleanField')(default=False),
                      keep_default=False)

        # Deleting field 'Sample.status'
        db.delete_column('labrack_sample', 'status')


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
            'Meta': {'object_name': 'Chassis'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentClassification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentClassification']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chassis_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.chemicalcomponent': {
            'Meta': {'object_name': 'ChemicalComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentClassification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentClassification']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chemicalcomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chemicalcomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chemicalcomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chemicalcomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.collection': {
            'Meta': {'object_name': 'Collection'},
            'chassis': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'dnaComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.NucleicAcidComponent']", 'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'proteinComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponent']", 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.componentclassification': {
            'Meta': {'object_name': 'ComponentClassification'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ComponentClassification']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200'})
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
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'})
        },
        'labrack.location': {
            'Meta': {'object_name': 'Location'},
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        'labrack.nucleicacidcomponent': {
            'Meta': {'object_name': 'NucleicAcidComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentClassification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentClassification']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.NucleicAcidComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.peptidecomponent': {
            'Meta': {'object_name': 'PeptideComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentClassification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentClassification']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'peptidecomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'peptidecomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'peptidecomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'peptidecomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.PeptideComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentClassification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentClassification']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'proteincomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'aliquotnb': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['labrack.Container']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sample_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'preparation_date': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime(2012, 6, 27, 0, 0)'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'ok'", 'max_length': '30'})
        },
        'labrack.samplecontent': {
            'Meta': {'object_name': 'SampleContent'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amount_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplecontent_amount_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentration_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplecontent_concentration_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sameplecontent'", 'to': "orm['labrack.Sample']"})
        },
        'labrack.samplepedigree': {
            'Meta': {'object_name': 'SamplePedigree'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amount_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplepedigree_amount_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentration_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplepedigree_concentration_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sameplepedigree'", 'to': "orm['labrack.Sample']"})
        },
        'labrack.selectivemarker': {
            'Meta': {'ordering': "('name',)", 'object_name': 'SelectiveMarker'},
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'})
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
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '10'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        }
    }

    complete_apps = ['labrack']