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
            ('createdBy', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='location_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('temperature', self.gf('django.db.models.fields.FloatField')()),
            ('room', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('registrationDate', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 6, 10, 0, 0))),
        ))
        db.send_create_signal('labrack', ['Location'])

        # Adding model 'Rack'
        db.create_table('labrack_rack', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('createdBy', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='rack_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('currentLocation', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Location'], null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('registrationDate', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 6, 10, 0, 0))),
        ))
        db.send_create_signal('labrack', ['Rack'])

        # Adding model 'Container'
        db.create_table('labrack_container', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('createdBy', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='container_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('containerType', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('rack', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Rack'])),
            ('registrationDate', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 6, 10, 0, 0))),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['Container'])

        # Adding model 'Unit'
        db.create_table('labrack_unit', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('conversion', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('unitType', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('labrack', ['Unit'])

        # Adding model 'Attachement'
        db.create_table('labrack_attachement', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uploadedFile', self.gf('django.db.models.fields.files.FileField')(max_length=100, null=True, blank=True)),
            ('fileType', self.gf('django.db.models.fields.CharField')(default='Type', max_length=30)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['Attachement'])

        # Adding model 'Source'
        db.create_table('labrack_source', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['Source'])

        # Adding model 'Component'
        db.create_table('labrack_component', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('createdBy', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='component_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('registrationDate', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 6, 10, 0, 0))),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('source', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Source'], null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
        ))
        db.send_create_signal('labrack', ['Component'])

        # Adding M2M table for field attachements on 'Component'
        db.create_table('labrack_component_attachements', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('attachement', models.ForeignKey(orm['labrack.attachement'], null=False))
        ))
        db.create_unique('labrack_component_attachements', ['component_id', 'attachement_id'])

        # Adding model 'CellType'
        db.create_table('labrack_celltype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('subTypeOf', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='subTypes', null=True, to=orm['labrack.CellType'])),
        ))
        db.send_create_signal('labrack', ['CellType'])

        # Adding model 'DnaComponentType'
        db.create_table('labrack_dnacomponenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('subTypeOf', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='subTypes', null=True, to=orm['labrack.DnaComponentType'])),
        ))
        db.send_create_signal('labrack', ['DnaComponentType'])

        # Adding model 'DnaComponent'
        db.create_table('labrack_dnacomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('circular', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('vectorBackbone', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.DnaComponent'], null=True, blank=True)),
            ('insert', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='Insert', null=True, to=orm['labrack.DnaComponent'])),
        ))
        db.send_create_signal('labrack', ['DnaComponent'])

        # Adding M2M table for field componentSubType on 'DnaComponent'
        db.create_table('labrack_dnacomponent_componentSubType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('dnacomponenttype', models.ForeignKey(orm['labrack.dnacomponenttype'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_componentSubType', ['dnacomponent_id', 'dnacomponenttype_id'])

        # Adding M2M table for field marker on 'DnaComponent'
        db.create_table('labrack_dnacomponent_marker', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('to_dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_marker', ['from_dnacomponent_id', 'to_dnacomponent_id'])

        # Adding model 'Sample'
        db.create_table('labrack_sample', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('createdBy', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='sample_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('container', self.gf('django.db.models.fields.related.ForeignKey')(related_name='samples', to=orm['labrack.Container'])),
            ('aliquotNr', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='ok', max_length=30)),
            ('comment', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('preparationDate', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 6, 10, 0, 0))),
            ('registrationDate', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 6, 10, 0, 0))),
            ('historyDescription', self.gf('django.db.models.fields.TextField')(max_length=255, null=True, blank=True)),
            ('source', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Source'], null=True, blank=True)),
            ('solvent', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
            ('concentration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('concentrationUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='concUnit', null=True, to=orm['labrack.Unit'])),
            ('amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amountUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='amouUnit', null=True, to=orm['labrack.Unit'])),
        ))
        db.send_create_signal('labrack', ['Sample'])

        # Adding M2M table for field attachements on 'Sample'
        db.create_table('labrack_sample_attachements', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sample', models.ForeignKey(orm['labrack.sample'], null=False)),
            ('attachement', models.ForeignKey(orm['labrack.attachement'], null=False))
        ))
        db.create_unique('labrack_sample_attachements', ['sample_id', 'attachement_id'])

        # Adding model 'Background'
        db.create_table('labrack_background', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=100)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['Background'])

        # Adding model 'Cell'
        db.create_table('labrack_cell', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('labrack', ['Cell'])

        # Adding M2M table for field componentSubType on 'Cell'
        db.create_table('labrack_cell_componentSubType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('cell', models.ForeignKey(orm['labrack.cell'], null=False)),
            ('celltype', models.ForeignKey(orm['labrack.celltype'], null=False))
        ))
        db.create_unique('labrack_cell_componentSubType', ['cell_id', 'celltype_id'])

        # Adding model 'PlasmidSample'
        db.create_table('labrack_plasmidsample', (
            ('sample_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Sample'], unique=True, primary_key=True)),
            ('dnaComponent', self.gf('django.db.models.fields.related.ForeignKey')(related_name='dna', to=orm['labrack.DnaComponent'])),
            ('inCell', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='in_Cell', null=True, to=orm['labrack.Cell'])),
        ))
        db.send_create_signal('labrack', ['PlasmidSample'])

        # Adding model 'CellSample'
        db.create_table('labrack_cellsample', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('cellComponent', self.gf('django.db.models.fields.related.ForeignKey')(related_name='cell', to=orm['labrack.Cell'])),
        ))
        db.send_create_signal('labrack', ['CellSample'])


    def backwards(self, orm):
        # Deleting model 'Location'
        db.delete_table('labrack_location')

        # Deleting model 'Rack'
        db.delete_table('labrack_rack')

        # Deleting model 'Container'
        db.delete_table('labrack_container')

        # Deleting model 'Unit'
        db.delete_table('labrack_unit')

        # Deleting model 'Attachement'
        db.delete_table('labrack_attachement')

        # Deleting model 'Source'
        db.delete_table('labrack_source')

        # Deleting model 'Component'
        db.delete_table('labrack_component')

        # Removing M2M table for field attachements on 'Component'
        db.delete_table('labrack_component_attachements')

        # Deleting model 'CellType'
        db.delete_table('labrack_celltype')

        # Deleting model 'DnaComponentType'
        db.delete_table('labrack_dnacomponenttype')

        # Deleting model 'DnaComponent'
        db.delete_table('labrack_dnacomponent')

        # Removing M2M table for field componentSubType on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_componentSubType')

        # Removing M2M table for field marker on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_marker')

        # Deleting model 'Sample'
        db.delete_table('labrack_sample')

        # Removing M2M table for field attachements on 'Sample'
        db.delete_table('labrack_sample_attachements')

        # Deleting model 'Background'
        db.delete_table('labrack_background')

        # Deleting model 'Cell'
        db.delete_table('labrack_cell')

        # Removing M2M table for field componentSubType on 'Cell'
        db.delete_table('labrack_cell_componentSubType')

        # Deleting model 'PlasmidSample'
        db.delete_table('labrack_plasmidsample')

        # Deleting model 'CellSample'
        db.delete_table('labrack_cellsample')


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
        'labrack.attachement': {
            'Meta': {'object_name': 'Attachement'},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'fileType': ('django.db.models.fields.CharField', [], {'default': "'Type'", 'max_length': '30'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uploadedFile': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'})
        },
        'labrack.background': {
            'Meta': {'object_name': 'Background'},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        'labrack.cell': {
            'Meta': {'object_name': 'Cell', '_ormbases': ['labrack.Component']},
            'componentSubType': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'Type'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.CellType']"}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.cellsample': {
            'Meta': {'object_name': 'CellSample', '_ormbases': ['labrack.Component']},
            'cellComponent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cell'", 'to': "orm['labrack.Cell']"}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.celltype': {
            'Meta': {'object_name': 'CellType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'to': "orm['labrack.CellType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.component': {
            'Meta': {'object_name': 'Component'},
            'attachements': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'Attachements'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.Attachement']"}),
            'comment': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'createdBy': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'component_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'registrationDate': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 6, 10, 0, 0)'}),
            'source': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Source']", 'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.container': {
            'Meta': {'ordering': "('displayId',)", 'object_name': 'Container'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'createdBy': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'container_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'rack': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Rack']"}),
            'registrationDate': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 6, 10, 0, 0)'})
        },
        'labrack.dnacomponent': {
            'Meta': {'object_name': 'DnaComponent', '_ormbases': ['labrack.Component']},
            'circular': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'componentSubType': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'Type'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.DnaComponentType']"}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'insert': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'Insert'", 'null': 'True', 'to': "orm['labrack.DnaComponent']"}),
            'marker': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'marker_rel_+'", 'null': 'True', 'to': "orm['labrack.DnaComponent']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'vectorBackbone': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.DnaComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.dnacomponenttype': {
            'Meta': {'object_name': 'DnaComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'to': "orm['labrack.DnaComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.location': {
            'Meta': {'object_name': 'Location'},
            'createdBy': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'location_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'registrationDate': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 6, 10, 0, 0)'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        'labrack.plasmidsample': {
            'Meta': {'object_name': 'PlasmidSample', '_ormbases': ['labrack.Sample']},
            'dnaComponent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'dna'", 'to': "orm['labrack.DnaComponent']"}),
            'inCell': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'in_Cell'", 'null': 'True', 'to': "orm['labrack.Cell']"}),
            'sample_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Sample']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.rack': {
            'Meta': {'object_name': 'Rack'},
            'createdBy': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'rack_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'currentLocation': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Location']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'registrationDate': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 6, 10, 0, 0)'})
        },
        'labrack.sample': {
            'Meta': {'object_name': 'Sample'},
            'aliquotNr': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amountUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'amouUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'attachements': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Attachement']", 'null': 'True', 'blank': 'True'}),
            'comment': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentrationUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'concUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['labrack.Container']"}),
            'createdBy': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sample_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'historyDescription': ('django.db.models.fields.TextField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'preparationDate': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 6, 10, 0, 0)'}),
            'registrationDate': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 6, 10, 0, 0)'}),
            'solvent': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'source': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Source']", 'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'ok'", 'max_length': '30'})
        },
        'labrack.source': {
            'Meta': {'object_name': 'Source'},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'})
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