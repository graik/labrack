# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Removing M2M table for field variantOf on 'Component'
        db.delete_table('labrack_component_variantOf')

        # Adding M2M table for field variantOf on 'ChemicalComponent'
        db.create_table(u'labrack_chemicalcomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('to_chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False))
        ))
        db.create_unique(u'labrack_chemicalcomponent_variantOf', ['from_chemicalcomponent_id', 'to_chemicalcomponent_id'])

        # Adding M2M table for field variantOf on 'Chassis'
        db.create_table(u'labrack_chassis_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('to_chassis', models.ForeignKey(orm['labrack.chassis'], null=False))
        ))
        db.create_unique(u'labrack_chassis_variantOf', ['from_chassis_id', 'to_chassis_id'])

        # Adding M2M table for field variantOf on 'ProteinComponent'
        db.create_table(u'labrack_proteincomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('to_proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False))
        ))
        db.create_unique(u'labrack_proteincomponent_variantOf', ['from_proteincomponent_id', 'to_proteincomponent_id'])

        # Adding M2M table for field variantOf on 'DnaComponent'
        db.create_table(u'labrack_dnacomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('to_dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False))
        ))
        db.create_unique(u'labrack_dnacomponent_variantOf', ['from_dnacomponent_id', 'to_dnacomponent_id'])

        # Adding M2M table for field variantOf on 'PeptideComponent'
        db.create_table(u'labrack_peptidecomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_peptidecomponent', models.ForeignKey(orm['labrack.peptidecomponent'], null=False)),
            ('to_peptidecomponent', models.ForeignKey(orm['labrack.peptidecomponent'], null=False))
        ))
        db.create_unique(u'labrack_peptidecomponent_variantOf', ['from_peptidecomponent_id', 'to_peptidecomponent_id'])


    def backwards(self, orm):
        # Adding M2M table for field variantOf on 'Component'
        db.create_table(u'labrack_component_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('to_component', models.ForeignKey(orm['labrack.component'], null=False))
        ))
        db.create_unique(u'labrack_component_variantOf', ['from_component_id', 'to_component_id'])

        # Removing M2M table for field variantOf on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_variantOf')

        # Removing M2M table for field variantOf on 'Chassis'
        db.delete_table('labrack_chassis_variantOf')

        # Removing M2M table for field variantOf on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_variantOf')

        # Removing M2M table for field variantOf on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_variantOf')

        # Removing M2M table for field variantOf on 'PeptideComponent'
        db.delete_table('labrack_peptidecomponent_variantOf')


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        'labrack.chassis': {
            'Meta': {'object_name': 'Chassis', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChassisComponentType']", 'null': 'True', 'blank': 'True'}),
            u'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.chassiscomponenttype': {
            'Meta': {'object_name': 'ChassisComponentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ChassisComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.chassissample': {
            'Meta': {'ordering': "('container', 'displayId')", 'object_name': 'ChassisSample', '_ormbases': ['labrack.Sample']},
            'chassis': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'Cell'", 'to': "orm['labrack.Chassis']"}),
            'derivedFrom': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'derivedFromSample'", 'null': 'True', 'to': "orm['labrack.ChassisSample']"}),
            'provenanceType': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '30', 'null': 'True', 'blank': 'True'}),
            u'sample_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Sample']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.chemicalcomponent': {
            'Meta': {'object_name': 'ChemicalComponent', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponentType']", 'null': 'True', 'blank': 'True'}),
            u'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.chemicalcomponenttype': {
            'Meta': {'object_name': 'ChemicalComponentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.collection': {
            'Meta': {'object_name': 'Collection'},
            'components': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Component']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.component': {
            'Meta': {'object_name': 'Component'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'component_created_by'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'component_owners'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 13, 0, 0)'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.container': {
            'Meta': {'object_name': 'Container'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'container_created_by'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_owners'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['auth.User']"}),
            'rack': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Rack']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 13, 0, 0)'})
        },
        'labrack.dnacomponent': {
            'GenBankfile': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'Meta': {'object_name': 'DnaComponent', '_ormbases': ['labrack.Component']},
            'circular': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.DnaComponentType']", 'null': 'True', 'blank': 'True'}),
            u'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.DnaComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.dnacomponenttype': {
            'Meta': {'object_name': 'DnaComponentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.DnaComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.dnasample': {
            'Meta': {'ordering': "('container', 'displayId')", 'object_name': 'DnaSample', '_ormbases': ['labrack.Sample']},
            'derivedFrom': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'derivedFromSample'", 'null': 'True', 'to': "orm['labrack.DnaSample']"}),
            'dnaConstruct': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'dnaSample'", 'null': 'True', 'to': "orm['labrack.DnaComponent']"}),
            'dna_sequenced': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '30', 'null': 'True', 'blank': 'True'}),
            'inChassis': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chassisSample'", 'null': 'True', 'to': "orm['labrack.Chassis']"}),
            'provenanceType': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '30', 'null': 'True', 'blank': 'True'}),
            u'sample_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Sample']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.dnasequenceannotation': {
            'Meta': {'object_name': 'DnaSequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'componentAnnotated': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'annotatedForComponent'", 'null': 'True', 'to': "orm['labrack.DnaComponent']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1', 'null': 'True', 'blank': 'True'}),
            'subComponent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'subComponentOf'", 'to': "orm['labrack.DnaComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.document': {
            'Meta': {'object_name': 'Document'},
            'docfile': ('django.db.models.fields.files.FileField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        'labrack.location': {
            'Meta': {'object_name': 'Location'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'location_created_by'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'location_owners'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 13, 0, 0)'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        'labrack.peptidecomponent': {
            'Meta': {'object_name': 'PeptideComponent', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.PeptideComponentType']", 'null': 'True', 'blank': 'True'}),
            u'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.PeptideComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.peptidecomponenttype': {
            'Meta': {'object_name': 'PeptideComponentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.PeptideComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponent': {
            'GenBankfile': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'Meta': {'object_name': 'ProteinComponent', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponentType']", 'null': 'True', 'blank': 'True'}),
            u'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponenttype': {
            'Meta': {'object_name': 'ProteinComponentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ProteinComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteinsequenceannotation': {
            'Meta': {'object_name': 'ProteinSequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'componentAnnotated': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'annotatedForComponent'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'subComponent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'subComponentOf'", 'to': "orm['labrack.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.rack': {
            'Meta': {'object_name': 'Rack'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'rack_created_by'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'current_location': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Location']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'rack_owners'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 13, 0, 0)'})
        },
        'labrack.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'aliquotNr': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amountUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'amouUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentrationUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'concUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['labrack.Container']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sample_created_by'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'historyDescription': ('django.db.models.fields.TextField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_owners'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['auth.User']"}),
            'preparation_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 13, 0, 0)'}),
            'reference_status': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 13, 0, 0)'}),
            'sampleCollection': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SampleCollection']", 'null': 'True', 'blank': 'True'}),
            'solvent': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'ok'", 'max_length': '30'})
        },
        'labrack.samplecollection': {
            'Meta': {'object_name': 'SampleCollection'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplecollection_created_by'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '25'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'samplecollection_owners'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['auth.User']"})
        },
        'labrack.sampleprovenance': {
            'Meta': {'object_name': 'SampleProvenance'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'provenanceType': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sampleProvenance'", 'to': "orm['labrack.Sample']"}),
            'sample_source': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sampleSource'", 'null': 'True', 'to': "orm['labrack.Sample']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'})
        },
        'labrack.selectivemarker': {
            'Meta': {'ordering': "('name',)", 'object_name': 'SelectiveMarker'},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'})
        },
        'labrack.unit': {
            'Meta': {'object_name': 'Unit'},
            'conversion': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'unitType': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        }
    }

    complete_apps = ['labrack']