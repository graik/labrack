# encoding: utf-8
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models

class Migration(SchemaMigration):

    def forwards(self, orm):
        
        # Adding model 'StorageContainer'
        db.create_table('brickit_storagecontainer', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('containerType', self.gf('django.db.models.fields.CharField')(max_length=30)),
        ))
        db.send_create_signal('brickit', ['StorageContainer'])

        # Adding M2M table for field users on 'StorageContainer'
        db.create_table('brickit_storagecontainer_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('storagecontainer', models.ForeignKey(orm['brickit.storagecontainer'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('brickit_storagecontainer_users', ['storagecontainer_id', 'user_id'])

        # Adding model 'Sample'
        db.create_table('brickit_sample', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('container', self.gf('django.db.models.fields.related.ForeignKey')(related_name='samples', to=orm['brickit.StorageContainer'])),
            ('wellType', self.gf('django.db.models.fields.CharField')(default='tube', max_length=30)),
            ('sampleType', self.gf('django.db.models.fields.CharField')(default='dna', max_length=100)),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('created', self.gf('django.db.models.fields.DateField')(auto_now_add=True, blank=True)),
            ('concentration', self.gf('django.db.models.fields.DecimalField')(default=50, max_digits=6, decimal_places=2)),
            ('concentrationUnit', self.gf('django.db.models.fields.CharField')(default='mg/l', max_length=10)),
        ))
        db.send_create_signal('brickit', ['Sample'])

        # Adding unique constraint on 'Sample', fields ['displayId', 'container']
        db.create_unique('brickit_sample', ['displayId', 'container_id'])

        # Adding M2M table for field users on 'Sample'
        db.create_table('brickit_sample_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sample', models.ForeignKey(orm['brickit.sample'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('brickit_sample_users', ['sample_id', 'user_id'])

        # Adding model 'ComponentType'
        db.create_table('brickit_componenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(unique=True, max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal('brickit', ['ComponentType'])

        # Adding M2M table for field subTypeOf on 'ComponentType'
        db.create_table('brickit_componenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_componenttype', models.ForeignKey(orm['brickit.componenttype'], null=False)),
            ('to_componenttype', models.ForeignKey(orm['brickit.componenttype'], null=False))
        ))
        db.create_unique('brickit_componenttype_subTypeOf', ['from_componenttype_id', 'to_componenttype_id'])

        # Adding model 'DnaComponent'
        db.create_table('brickit_dnacomponent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('translatesTo', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='encodedBy', null=True, to=orm['brickit.ProteinComponent'])),
            ('optimizedFor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['brickit.Chassis'], null=True, blank=True)),
        ))
        db.send_create_signal('brickit', ['DnaComponent'])

        # Adding M2M table for field componentType on 'DnaComponent'
        db.create_table('brickit_dnacomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['brickit.dnacomponent'], null=False)),
            ('componenttype', models.ForeignKey(orm['brickit.componenttype'], null=False))
        ))
        db.create_unique('brickit_dnacomponent_componentType', ['dnacomponent_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'DnaComponent'
        db.create_table('brickit_dnacomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_dnacomponent', models.ForeignKey(orm['brickit.dnacomponent'], null=False)),
            ('to_dnacomponent', models.ForeignKey(orm['brickit.dnacomponent'], null=False))
        ))
        db.create_unique('brickit_dnacomponent_variantOf', ['from_dnacomponent_id', 'to_dnacomponent_id'])

        # Adding M2M table for field users on 'DnaComponent'
        db.create_table('brickit_dnacomponent_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['brickit.dnacomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('brickit_dnacomponent_users', ['dnacomponent_id', 'user_id'])

        # Adding model 'ProteinComponent'
        db.create_table('brickit_proteincomponent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal('brickit', ['ProteinComponent'])

        # Adding M2M table for field componentType on 'ProteinComponent'
        db.create_table('brickit_proteincomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['brickit.proteincomponent'], null=False)),
            ('componenttype', models.ForeignKey(orm['brickit.componenttype'], null=False))
        ))
        db.create_unique('brickit_proteincomponent_componentType', ['proteincomponent_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'ProteinComponent'
        db.create_table('brickit_proteincomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_proteincomponent', models.ForeignKey(orm['brickit.proteincomponent'], null=False)),
            ('to_proteincomponent', models.ForeignKey(orm['brickit.proteincomponent'], null=False))
        ))
        db.create_unique('brickit_proteincomponent_variantOf', ['from_proteincomponent_id', 'to_proteincomponent_id'])

        # Adding M2M table for field users on 'ProteinComponent'
        db.create_table('brickit_proteincomponent_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['brickit.proteincomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('brickit_proteincomponent_users', ['proteincomponent_id', 'user_id'])

        # Adding model 'Chassis'
        db.create_table('brickit_chassis', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal('brickit', ['Chassis'])

        # Adding M2M table for field componentType on 'Chassis'
        db.create_table('brickit_chassis_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['brickit.chassis'], null=False)),
            ('componenttype', models.ForeignKey(orm['brickit.componenttype'], null=False))
        ))
        db.create_unique('brickit_chassis_componentType', ['chassis_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'Chassis'
        db.create_table('brickit_chassis_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chassis', models.ForeignKey(orm['brickit.chassis'], null=False)),
            ('to_chassis', models.ForeignKey(orm['brickit.chassis'], null=False))
        ))
        db.create_unique('brickit_chassis_variantOf', ['from_chassis_id', 'to_chassis_id'])

        # Adding M2M table for field users on 'Chassis'
        db.create_table('brickit_chassis_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['brickit.chassis'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('brickit_chassis_users', ['chassis_id', 'user_id'])


    def backwards(self, orm):
        
        # Removing unique constraint on 'Sample', fields ['displayId', 'container']
        db.delete_unique('brickit_sample', ['displayId', 'container_id'])

        # Deleting model 'StorageContainer'
        db.delete_table('brickit_storagecontainer')

        # Removing M2M table for field users on 'StorageContainer'
        db.delete_table('brickit_storagecontainer_users')

        # Deleting model 'Sample'
        db.delete_table('brickit_sample')

        # Removing M2M table for field users on 'Sample'
        db.delete_table('brickit_sample_users')

        # Deleting model 'ComponentType'
        db.delete_table('brickit_componenttype')

        # Removing M2M table for field subTypeOf on 'ComponentType'
        db.delete_table('brickit_componenttype_subTypeOf')

        # Deleting model 'DnaComponent'
        db.delete_table('brickit_dnacomponent')

        # Removing M2M table for field componentType on 'DnaComponent'
        db.delete_table('brickit_dnacomponent_componentType')

        # Removing M2M table for field variantOf on 'DnaComponent'
        db.delete_table('brickit_dnacomponent_variantOf')

        # Removing M2M table for field users on 'DnaComponent'
        db.delete_table('brickit_dnacomponent_users')

        # Deleting model 'ProteinComponent'
        db.delete_table('brickit_proteincomponent')

        # Removing M2M table for field componentType on 'ProteinComponent'
        db.delete_table('brickit_proteincomponent_componentType')

        # Removing M2M table for field variantOf on 'ProteinComponent'
        db.delete_table('brickit_proteincomponent_variantOf')

        # Removing M2M table for field users on 'ProteinComponent'
        db.delete_table('brickit_proteincomponent_users')

        # Deleting model 'Chassis'
        db.delete_table('brickit_chassis')

        # Removing M2M table for field componentType on 'Chassis'
        db.delete_table('brickit_chassis_componentType')

        # Removing M2M table for field variantOf on 'Chassis'
        db.delete_table('brickit_chassis_variantOf')

        # Removing M2M table for field users on 'Chassis'
        db.delete_table('brickit_chassis_users')


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
        'brickit.chassis': {
            'Meta': {'object_name': 'Chassis'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.Chassis']", 'null': 'True', 'blank': 'True'})
        },
        'brickit.componenttype': {
            'Meta': {'object_name': 'ComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['brickit.ComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'unique': 'True', 'max_length': '200'})
        },
        'brickit.dnacomponent': {
            'Meta': {'object_name': 'DnaComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['brickit.Chassis']", 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['brickit.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.DnaComponent']", 'null': 'True', 'blank': 'True'})
        },
        'brickit.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ProteinComponent']", 'null': 'True', 'blank': 'True'})
        },
        'brickit.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '50', 'max_digits': '6', 'decimal_places': '2'}),
            'concentrationUnit': ('django.db.models.fields.CharField', [], {'default': "'mg/l'", 'max_length': '10'}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['brickit.StorageContainer']"}),
            'created': ('django.db.models.fields.DateField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sampleType': ('django.db.models.fields.CharField', [], {'default': "'dna'", 'max_length': '100'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'wellType': ('django.db.models.fields.CharField', [], {'default': "'tube'", 'max_length': '30'})
        },
        'brickit.storagecontainer': {
            'Meta': {'ordering': "('displayId',)", 'object_name': 'StorageContainer'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['auth.User']", 'null': 'True', 'db_index': 'True', 'blank': 'True'})
        },
        'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['brickit']
