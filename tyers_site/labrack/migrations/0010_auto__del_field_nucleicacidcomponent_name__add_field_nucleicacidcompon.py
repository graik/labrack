# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting field 'NucleicAcidComponent.name'
        db.delete_column('labrack_nucleicacidcomponent', 'name')

        # Adding field 'NucleicAcidComponent.created_by'
        db.add_column('labrack_nucleicacidcomponent', 'created_by',
                      self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='nucleicacidcomponent_created_by', null=True, to=orm['auth.User']),
                      keep_default=False)

        # Removing M2M table for field users on 'NucleicAcidComponent'
        db.delete_table('labrack_nucleicacidcomponent_users')

        # Adding M2M table for field owners on 'NucleicAcidComponent'
        db.create_table('labrack_nucleicacidcomponent_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('nucleicacidcomponent', models.ForeignKey(orm['labrack.nucleicacidcomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_nucleicacidcomponent_owners', ['nucleicacidcomponent_id', 'user_id'])

        # Adding M2M table for field group_read on 'NucleicAcidComponent'
        db.create_table('labrack_nucleicacidcomponent_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('nucleicacidcomponent', models.ForeignKey(orm['labrack.nucleicacidcomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_nucleicacidcomponent_group_read', ['nucleicacidcomponent_id', 'group_id'])

        # Adding M2M table for field group_write on 'NucleicAcidComponent'
        db.create_table('labrack_nucleicacidcomponent_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('nucleicacidcomponent', models.ForeignKey(orm['labrack.nucleicacidcomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_nucleicacidcomponent_group_write', ['nucleicacidcomponent_id', 'group_id'])


        # Changing field 'NucleicAcidComponent.shortDescription'
        db.alter_column('labrack_nucleicacidcomponent', 'shortDescription', self.gf('django.db.models.fields.CharField')(max_length=50))
        # Deleting field 'Chassis.name'
        db.delete_column('labrack_chassis', 'name')

        # Adding field 'Chassis.created_by'
        db.add_column('labrack_chassis', 'created_by',
                      self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='chassis_created_by', null=True, to=orm['auth.User']),
                      keep_default=False)

        # Removing M2M table for field users on 'Chassis'
        db.delete_table('labrack_chassis_users')

        # Adding M2M table for field owners on 'Chassis'
        db.create_table('labrack_chassis_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_chassis_owners', ['chassis_id', 'user_id'])

        # Adding M2M table for field group_read on 'Chassis'
        db.create_table('labrack_chassis_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_chassis_group_read', ['chassis_id', 'group_id'])

        # Adding M2M table for field group_write on 'Chassis'
        db.create_table('labrack_chassis_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_chassis_group_write', ['chassis_id', 'group_id'])


        # Changing field 'Chassis.shortDescription'
        db.alter_column('labrack_chassis', 'shortDescription', self.gf('django.db.models.fields.CharField')(max_length=50))
        # Deleting field 'ProteinComponent.name'
        db.delete_column('labrack_proteincomponent', 'name')

        # Adding field 'ProteinComponent.created_by'
        db.add_column('labrack_proteincomponent', 'created_by',
                      self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='proteincomponent_created_by', null=True, to=orm['auth.User']),
                      keep_default=False)

        # Removing M2M table for field users on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_users')

        # Adding M2M table for field owners on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_owners', ['proteincomponent_id', 'user_id'])

        # Adding M2M table for field group_read on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_group_read', ['proteincomponent_id', 'group_id'])

        # Adding M2M table for field group_write on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_group_write', ['proteincomponent_id', 'group_id'])


        # Changing field 'ProteinComponent.shortDescription'
        db.alter_column('labrack_proteincomponent', 'shortDescription', self.gf('django.db.models.fields.CharField')(max_length=50))

    def backwards(self, orm):
        # Adding field 'NucleicAcidComponent.name'
        db.add_column('labrack_nucleicacidcomponent', 'name',
                      self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True),
                      keep_default=False)

        # Deleting field 'NucleicAcidComponent.created_by'
        db.delete_column('labrack_nucleicacidcomponent', 'created_by_id')

        # Adding M2M table for field users on 'NucleicAcidComponent'
        db.create_table('labrack_nucleicacidcomponent_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('nucleicacidcomponent', models.ForeignKey(orm['labrack.nucleicacidcomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_nucleicacidcomponent_users', ['nucleicacidcomponent_id', 'user_id'])

        # Removing M2M table for field owners on 'NucleicAcidComponent'
        db.delete_table('labrack_nucleicacidcomponent_owners')

        # Removing M2M table for field group_read on 'NucleicAcidComponent'
        db.delete_table('labrack_nucleicacidcomponent_group_read')

        # Removing M2M table for field group_write on 'NucleicAcidComponent'
        db.delete_table('labrack_nucleicacidcomponent_group_write')


        # Changing field 'NucleicAcidComponent.shortDescription'
        db.alter_column('labrack_nucleicacidcomponent', 'shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200))
        # Adding field 'Chassis.name'
        db.add_column('labrack_chassis', 'name',
                      self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True),
                      keep_default=False)

        # Deleting field 'Chassis.created_by'
        db.delete_column('labrack_chassis', 'created_by_id')

        # Adding M2M table for field users on 'Chassis'
        db.create_table('labrack_chassis_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_chassis_users', ['chassis_id', 'user_id'])

        # Removing M2M table for field owners on 'Chassis'
        db.delete_table('labrack_chassis_owners')

        # Removing M2M table for field group_read on 'Chassis'
        db.delete_table('labrack_chassis_group_read')

        # Removing M2M table for field group_write on 'Chassis'
        db.delete_table('labrack_chassis_group_write')


        # Changing field 'Chassis.shortDescription'
        db.alter_column('labrack_chassis', 'shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200))
        # Adding field 'ProteinComponent.name'
        db.add_column('labrack_proteincomponent', 'name',
                      self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True),
                      keep_default=False)

        # Deleting field 'ProteinComponent.created_by'
        db.delete_column('labrack_proteincomponent', 'created_by_id')

        # Adding M2M table for field users on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_users', ['proteincomponent_id', 'user_id'])

        # Removing M2M table for field owners on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_owners')

        # Removing M2M table for field group_read on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_group_read')

        # Removing M2M table for field group_write on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_group_write')


        # Changing field 'ProteinComponent.shortDescription'
        db.alter_column('labrack_proteincomponent', 'shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200))

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
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chassis_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'})
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
        'labrack.componenttype': {
            'Meta': {'object_name': 'ComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ComponentType']"}),
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
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'nucleicacidcomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.NucleicAcidComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'proteincomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
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
            'preparation_date': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime(2012, 6, 21, 0, 0)'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'})
        },
        'labrack.samplecontent': {
            'Meta': {'object_name': 'SampleContent'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amount_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'amount_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentration_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'concentration_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Sample']"})
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
            'name': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        }
    }

    complete_apps = ['labrack']