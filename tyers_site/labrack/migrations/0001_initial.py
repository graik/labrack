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
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('temperature', self.gf('django.db.models.fields.FloatField')()),
            ('room', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['Location'])

        # Adding model 'Container'
        db.create_table('labrack_container', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('containerType', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('location', self.gf('django.db.models.fields.related.ForeignKey')(related_name='containers', to=orm['labrack.Location'])),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='container_created_by', null=True, to=orm['auth.User'])),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['Container'])

        # Adding M2M table for field owners on 'Container'
        db.create_table('labrack_container_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['labrack.container'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_container_owners', ['container_id', 'user_id'])

        # Adding M2M table for field group_read on 'Container'
        db.create_table('labrack_container_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['labrack.container'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_container_group_read', ['container_id', 'group_id'])

        # Adding M2M table for field group_write on 'Container'
        db.create_table('labrack_container_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('container', models.ForeignKey(orm['labrack.container'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_container_group_write', ['container_id', 'group_id'])

        # Adding model 'Unit'
        db.create_table('labrack_unit', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('conversion', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('unitType', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('labrack', ['Unit'])

        # Adding model 'Sample'
        db.create_table('labrack_sample', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('container', self.gf('django.db.models.fields.related.ForeignKey')(related_name='samples', to=orm['labrack.Container'])),
            ('aliquotNr', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='ok', max_length=30)),
            ('attachment', self.gf('django.db.models.fields.files.FileField')(max_length=100, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('preparation_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2012, 7, 13, 0, 0))),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='sample_created_by', null=True, to=orm['auth.User'])),
        ))
        db.send_create_signal('labrack', ['Sample'])

        # Adding unique constraint on 'Sample', fields ['displayId', 'container']
        db.create_unique('labrack_sample', ['displayId', 'container_id'])

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

        # Adding model 'Project'
        db.create_table('labrack_project', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=25)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='project_created_by', null=True, to=orm['auth.User'])),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['Project'])

        # Adding M2M table for field owners on 'Project'
        db.create_table('labrack_project_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('project', models.ForeignKey(orm['labrack.project'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_project_owners', ['project_id', 'user_id'])

        # Adding M2M table for field group_read on 'Project'
        db.create_table('labrack_project_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('project', models.ForeignKey(orm['labrack.project'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_project_group_read', ['project_id', 'group_id'])

        # Adding M2M table for field group_write on 'Project'
        db.create_table('labrack_project_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('project', models.ForeignKey(orm['labrack.project'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_project_group_write', ['project_id', 'group_id'])

        # Adding model 'SampleLink'
        db.create_table('labrack_samplelink', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample', self.gf('django.db.models.fields.related.ForeignKey')(related_name='sameplelink', to=orm['labrack.Sample'])),
            ('linkType', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('link', self.gf('django.db.models.fields.CharField')(max_length=1000)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['SampleLink'])

        # Adding model 'SampleContent'
        db.create_table('labrack_samplecontent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample', self.gf('django.db.models.fields.related.ForeignKey')(related_name='samepleContent', to=orm['labrack.Sample'])),
            ('content_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['contenttypes.ContentType'])),
            ('object_id', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amountUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='amountUnit', null=True, to=orm['labrack.Unit'])),
            ('concentration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('concentrationUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='concentrationUnit', null=True, to=orm['labrack.Unit'])),
        ))
        db.send_create_signal('labrack', ['SampleContent'])

        # Adding model 'SamplePedigree'
        db.create_table('labrack_samplepedigree', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample', self.gf('django.db.models.fields.related.ForeignKey')(related_name='sameplepedigree', to=orm['labrack.Sample'])),
            ('sample_source', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Sample'], null=True)),
            ('amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amountUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='samplepedigree_amount_unit', null=True, to=orm['labrack.Unit'])),
            ('concentration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('concentrationUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='samplepedigree_concentration_unit', null=True, to=orm['labrack.Unit'])),
        ))
        db.send_create_signal('labrack', ['SamplePedigree'])

        # Adding model 'ComponentType'
        db.create_table('labrack_componenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('labrack', ['ComponentType'])

        # Adding M2M table for field subTypeOf on 'ComponentType'
        db.create_table('labrack_componenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False)),
            ('to_componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_componenttype_subTypeOf', ['from_componenttype_id', 'to_componenttype_id'])

        # Adding model 'Component'
        db.create_table('labrack_component', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, null=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='component_created_by', null=True, to=orm['auth.User'])),
        ))
        db.send_create_signal('labrack', ['Component'])

        # Adding M2M table for field componentType on 'Component'
        db.create_table('labrack_component_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_component_componentType', ['component_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'Component'
        db.create_table('labrack_component_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('to_component', models.ForeignKey(orm['labrack.component'], null=False))
        ))
        db.create_unique('labrack_component_variantOf', ['from_component_id', 'to_component_id'])

        # Adding M2M table for field annotations on 'Component'
        db.create_table('labrack_component_annotations', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False))
        ))
        db.create_unique('labrack_component_annotations', ['component_id', 'sequenceannotation_id'])

        # Adding M2M table for field owners on 'Component'
        db.create_table('labrack_component_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_component_owners', ['component_id', 'user_id'])

        # Adding M2M table for field group_read on 'Component'
        db.create_table('labrack_component_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_component_group_read', ['component_id', 'group_id'])

        # Adding M2M table for field group_write on 'Component'
        db.create_table('labrack_component_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_component_group_write', ['component_id', 'group_id'])

        # Adding model 'DnaComponent'
        db.create_table('labrack_dnacomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('translatesTo', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='encodedBy', null=True, to=orm['labrack.ProteinComponent'])),
            ('optimizedFor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Chassis'], null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['DnaComponent'])

        # Adding model 'Chassis'
        db.create_table('labrack_chassis', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('labrack', ['Chassis'])

        # Adding model 'ProteinComponent'
        db.create_table('labrack_proteincomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['ProteinComponent'])

        # Adding model 'PeptideComponent'
        db.create_table('labrack_peptidecomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['PeptideComponent'])

        # Adding model 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('labrack', ['ChemicalComponent'])

        # Adding model 'SelectiveMarker'
        db.create_table('labrack_selectivemarker', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['SelectiveMarker'])

        # Adding model 'SequenceAnnotation'
        db.create_table('labrack_sequenceannotation', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('bioStart', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('bioEnd', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('strand', self.gf('django.db.models.fields.CharField')(max_length=1, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['SequenceAnnotation'])

        # Adding M2M table for field precedes on 'SequenceAnnotation'
        db.create_table('labrack_sequenceannotation_precedes', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False)),
            ('to_sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False))
        ))
        db.create_unique('labrack_sequenceannotation_precedes', ['from_sequenceannotation_id', 'to_sequenceannotation_id'])

        # Adding model 'Collection'
        db.create_table('labrack_collection', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['Collection'])

        # Adding M2M table for field components on 'Collection'
        db.create_table('labrack_collection_components', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('collection', models.ForeignKey(orm['labrack.collection'], null=False)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False))
        ))
        db.create_unique('labrack_collection_components', ['collection_id', 'component_id'])


    def backwards(self, orm):
        # Removing unique constraint on 'Sample', fields ['displayId', 'container']
        db.delete_unique('labrack_sample', ['displayId', 'container_id'])

        # Deleting model 'Location'
        db.delete_table('labrack_location')

        # Deleting model 'Container'
        db.delete_table('labrack_container')

        # Removing M2M table for field owners on 'Container'
        db.delete_table('labrack_container_owners')

        # Removing M2M table for field group_read on 'Container'
        db.delete_table('labrack_container_group_read')

        # Removing M2M table for field group_write on 'Container'
        db.delete_table('labrack_container_group_write')

        # Deleting model 'Unit'
        db.delete_table('labrack_unit')

        # Deleting model 'Sample'
        db.delete_table('labrack_sample')

        # Removing M2M table for field owners on 'Sample'
        db.delete_table('labrack_sample_owners')

        # Removing M2M table for field group_read on 'Sample'
        db.delete_table('labrack_sample_group_read')

        # Removing M2M table for field group_write on 'Sample'
        db.delete_table('labrack_sample_group_write')

        # Deleting model 'Project'
        db.delete_table('labrack_project')

        # Removing M2M table for field owners on 'Project'
        db.delete_table('labrack_project_owners')

        # Removing M2M table for field group_read on 'Project'
        db.delete_table('labrack_project_group_read')

        # Removing M2M table for field group_write on 'Project'
        db.delete_table('labrack_project_group_write')

        # Deleting model 'SampleLink'
        db.delete_table('labrack_samplelink')

        # Deleting model 'SampleContent'
        db.delete_table('labrack_samplecontent')

        # Deleting model 'SamplePedigree'
        db.delete_table('labrack_samplepedigree')

        # Deleting model 'ComponentType'
        db.delete_table('labrack_componenttype')

        # Removing M2M table for field subTypeOf on 'ComponentType'
        db.delete_table('labrack_componenttype_subTypeOf')

        # Deleting model 'Component'
        db.delete_table('labrack_component')

        # Removing M2M table for field componentType on 'Component'
        db.delete_table('labrack_component_componentType')

        # Removing M2M table for field variantOf on 'Component'
        db.delete_table('labrack_component_variantOf')

        # Removing M2M table for field annotations on 'Component'
        db.delete_table('labrack_component_annotations')

        # Removing M2M table for field owners on 'Component'
        db.delete_table('labrack_component_owners')

        # Removing M2M table for field group_read on 'Component'
        db.delete_table('labrack_component_group_read')

        # Removing M2M table for field group_write on 'Component'
        db.delete_table('labrack_component_group_write')

        # Deleting model 'DnaComponent'
        db.delete_table('labrack_dnacomponent')

        # Deleting model 'Chassis'
        db.delete_table('labrack_chassis')

        # Deleting model 'ProteinComponent'
        db.delete_table('labrack_proteincomponent')

        # Deleting model 'PeptideComponent'
        db.delete_table('labrack_peptidecomponent')

        # Deleting model 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent')

        # Deleting model 'SelectiveMarker'
        db.delete_table('labrack_selectivemarker')

        # Deleting model 'SequenceAnnotation'
        db.delete_table('labrack_sequenceannotation')

        # Removing M2M table for field precedes on 'SequenceAnnotation'
        db.delete_table('labrack_sequenceannotation_precedes')

        # Deleting model 'Collection'
        db.delete_table('labrack_collection')

        # Removing M2M table for field components on 'Collection'
        db.delete_table('labrack_collection_components')


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
        'labrack.project': {
            'Meta': {'object_name': 'Project'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'project_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'project_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'project_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '25'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'project_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'})
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
            'preparation_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2012, 7, 13, 0, 0)'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'ok'", 'max_length': '30'})
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
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samepleContent'", 'to': "orm['labrack.Sample']"})
        },
        'labrack.samplelink': {
            'Meta': {'object_name': 'SampleLink'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'link': ('django.db.models.fields.CharField', [], {'max_length': '1000'}),
            'linkType': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sameplelink'", 'to': "orm['labrack.Sample']"}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'})
        },
        'labrack.samplepedigree': {
            'Meta': {'object_name': 'SamplePedigree'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amountUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplepedigree_amount_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentrationUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplepedigree_concentration_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sameplepedigree'", 'to': "orm['labrack.Sample']"}),
            'sample_source': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Sample']", 'null': 'True'})
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