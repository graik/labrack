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
            ('empty', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('preparation_date', self.gf('django.db.models.fields.DateTimeField')(default=datetime.datetime(2012, 7, 5, 0, 0))),
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

        # Adding model 'SampleContent'
        db.create_table('labrack_samplecontent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample', self.gf('django.db.models.fields.related.ForeignKey')(related_name='sameplecontent', to=orm['labrack.Sample'])),
            ('content_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['contenttypes.ContentType'])),
            ('object_id', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amount_unit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='mass_unit', null=True, to=orm['labrack.Unit'])),
            ('concentration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('concentration_unit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='concentration_unit', null=True, to=orm['labrack.Unit'])),
        ))
        db.send_create_signal('labrack', ['SampleContent'])

        # Adding model 'ComponentType'
        db.create_table('labrack_componenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
        ))
        db.send_create_signal('labrack', ['ComponentType'])

        # Adding M2M table for field subTypeOf on 'ComponentType'
        db.create_table('labrack_componenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False)),
            ('to_componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_componenttype_subTypeOf', ['from_componenttype_id', 'to_componenttype_id'])

        # Adding model 'DnaComponent'
        db.create_table('labrack_dnacomponent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, null=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='dnacomponent_created_by', null=True, to=orm['auth.User'])),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('translatesTo', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='encodedBy', null=True, to=orm['labrack.ProteinComponent'])),
            ('optimizedFor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Chassis'], null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['DnaComponent'])

        # Adding M2M table for field componentType on 'DnaComponent'
        db.create_table('labrack_dnacomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_componentType', ['dnacomponent_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'DnaComponent'
        db.create_table('labrack_dnacomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('to_dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_variantOf', ['from_dnacomponent_id', 'to_dnacomponent_id'])

        # Adding M2M table for field annotations on 'DnaComponent'
        db.create_table('labrack_dnacomponent_annotations', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_annotations', ['dnacomponent_id', 'sequenceannotation_id'])

        # Adding M2M table for field owners on 'DnaComponent'
        db.create_table('labrack_dnacomponent_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_owners', ['dnacomponent_id', 'user_id'])

        # Adding M2M table for field group_read on 'DnaComponent'
        db.create_table('labrack_dnacomponent_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_group_read', ['dnacomponent_id', 'group_id'])

        # Adding M2M table for field group_write on 'DnaComponent'
        db.create_table('labrack_dnacomponent_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_group_write', ['dnacomponent_id', 'group_id'])

        # Adding model 'Chassis'
        db.create_table('labrack_chassis', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, null=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='chassis_created_by', null=True, to=orm['auth.User'])),
        ))
        db.send_create_signal('labrack', ['Chassis'])

        # Adding M2M table for field componentType on 'Chassis'
        db.create_table('labrack_chassis_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_chassis_componentType', ['chassis_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'Chassis'
        db.create_table('labrack_chassis_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('to_chassis', models.ForeignKey(orm['labrack.chassis'], null=False))
        ))
        db.create_unique('labrack_chassis_variantOf', ['from_chassis_id', 'to_chassis_id'])

        # Adding M2M table for field annotations on 'Chassis'
        db.create_table('labrack_chassis_annotations', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False))
        ))
        db.create_unique('labrack_chassis_annotations', ['chassis_id', 'sequenceannotation_id'])

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

        # Adding model 'ProteinComponent'
        db.create_table('labrack_proteincomponent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, null=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='proteincomponent_created_by', null=True, to=orm['auth.User'])),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['ProteinComponent'])

        # Adding M2M table for field componentType on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_componentType', ['proteincomponent_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('to_proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_variantOf', ['from_proteincomponent_id', 'to_proteincomponent_id'])

        # Adding M2M table for field annotations on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_annotations', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_annotations', ['proteincomponent_id', 'sequenceannotation_id'])

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

        # Adding model 'PeptideComponent'
        db.create_table('labrack_peptidecomponent', (
            ('proteincomponent_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.ProteinComponent'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('labrack', ['PeptideComponent'])

        # Adding model 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
            ('abstract', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, null=True, blank=True)),
            ('modification_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='chemicalcomponent_created_by', null=True, to=orm['auth.User'])),
        ))
        db.send_create_signal('labrack', ['ChemicalComponent'])

        # Adding M2M table for field componentType on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('componenttype', models.ForeignKey(orm['labrack.componenttype'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_componentType', ['chemicalcomponent_id', 'componenttype_id'])

        # Adding M2M table for field variantOf on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('to_chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_variantOf', ['from_chemicalcomponent_id', 'to_chemicalcomponent_id'])

        # Adding M2M table for field annotations on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_annotations', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('sequenceannotation', models.ForeignKey(orm['labrack.sequenceannotation'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_annotations', ['chemicalcomponent_id', 'sequenceannotation_id'])

        # Adding M2M table for field owners on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_owners', ['chemicalcomponent_id', 'user_id'])

        # Adding M2M table for field group_read on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_group_read', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_group_read', ['chemicalcomponent_id', 'group_id'])

        # Adding M2M table for field group_write on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_group_write', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('group', models.ForeignKey(orm['auth.group'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_group_write', ['chemicalcomponent_id', 'group_id'])

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

        # Adding M2M table for field dnaComponents on 'Collection'
        db.create_table('labrack_collection_dnaComponents', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('collection', models.ForeignKey(orm['labrack.collection'], null=False)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False))
        ))
        db.create_unique('labrack_collection_dnaComponents', ['collection_id', 'dnacomponent_id'])

        # Adding M2M table for field proteinComponents on 'Collection'
        db.create_table('labrack_collection_proteinComponents', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('collection', models.ForeignKey(orm['labrack.collection'], null=False)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False))
        ))
        db.create_unique('labrack_collection_proteinComponents', ['collection_id', 'proteincomponent_id'])

        # Adding M2M table for field chassis on 'Collection'
        db.create_table('labrack_collection_chassis', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('collection', models.ForeignKey(orm['labrack.collection'], null=False)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False))
        ))
        db.create_unique('labrack_collection_chassis', ['collection_id', 'chassis_id'])


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

        # Deleting model 'SampleContent'
        db.delete_table('labrack_samplecontent')

        # Deleting model 'ComponentType'
        db.delete_table('labrack_componenttype')

        # Removing M2M table for field subTypeOf on 'ComponentType'
        db.delete_table('labrack_componenttype_subTypeOf')

        # Deleting model 'DnaComponent'
        db.delete_table('labrack_dnacomponent')

        # Removing M2M table for field componentType on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_componentType')

        # Removing M2M table for field variantOf on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_variantOf')

        # Removing M2M table for field annotations on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_annotations')

        # Removing M2M table for field owners on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_owners')

        # Removing M2M table for field group_read on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_group_read')

        # Removing M2M table for field group_write on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_group_write')

        # Deleting model 'Chassis'
        db.delete_table('labrack_chassis')

        # Removing M2M table for field componentType on 'Chassis'
        db.delete_table('labrack_chassis_componentType')

        # Removing M2M table for field variantOf on 'Chassis'
        db.delete_table('labrack_chassis_variantOf')

        # Removing M2M table for field annotations on 'Chassis'
        db.delete_table('labrack_chassis_annotations')

        # Removing M2M table for field owners on 'Chassis'
        db.delete_table('labrack_chassis_owners')

        # Removing M2M table for field group_read on 'Chassis'
        db.delete_table('labrack_chassis_group_read')

        # Removing M2M table for field group_write on 'Chassis'
        db.delete_table('labrack_chassis_group_write')

        # Deleting model 'ProteinComponent'
        db.delete_table('labrack_proteincomponent')

        # Removing M2M table for field componentType on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_componentType')

        # Removing M2M table for field variantOf on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_variantOf')

        # Removing M2M table for field annotations on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_annotations')

        # Removing M2M table for field owners on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_owners')

        # Removing M2M table for field group_read on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_group_read')

        # Removing M2M table for field group_write on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_group_write')

        # Deleting model 'PeptideComponent'
        db.delete_table('labrack_peptidecomponent')

        # Deleting model 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent')

        # Removing M2M table for field componentType on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_componentType')

        # Removing M2M table for field variantOf on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_variantOf')

        # Removing M2M table for field annotations on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_annotations')

        # Removing M2M table for field owners on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_owners')

        # Removing M2M table for field group_read on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_group_read')

        # Removing M2M table for field group_write on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_group_write')

        # Deleting model 'SelectiveMarker'
        db.delete_table('labrack_selectivemarker')

        # Deleting model 'SequenceAnnotation'
        db.delete_table('labrack_sequenceannotation')

        # Removing M2M table for field precedes on 'SequenceAnnotation'
        db.delete_table('labrack_sequenceannotation_precedes')

        # Deleting model 'Collection'
        db.delete_table('labrack_collection')

        # Removing M2M table for field dnaComponents on 'Collection'
        db.delete_table('labrack_collection_dnaComponents')

        # Removing M2M table for field proteinComponents on 'Collection'
        db.delete_table('labrack_collection_proteinComponents')

        # Removing M2M table for field chassis on 'Collection'
        db.delete_table('labrack_collection_chassis')


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
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chassis_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.chemicalcomponent': {
            'Meta': {'object_name': 'ChemicalComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chemicalcomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chemicalcomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chemicalcomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'chemicalcomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.collection': {
            'Meta': {'object_name': 'Collection'},
            'chassis': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'dnaComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.DnaComponent']", 'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'proteinComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponent']", 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.componenttype': {
            'Meta': {'object_name': 'ComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
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
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"})
        },
        'labrack.dnacomponent': {
            'Meta': {'object_name': 'DnaComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'dnacomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'dnacomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'dnacomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'dnacomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.DnaComponent']", 'null': 'True', 'blank': 'True'})
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
            'Meta': {'object_name': 'PeptideComponent', '_ormbases': ['labrack.ProteinComponent']},
            'proteincomponent_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.ProteinComponent']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'proteincomponent_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'proteincomponent_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponent']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'aliquotNr': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['labrack.Container']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sample_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'empty': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'group_read': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_groups_read'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'group_write': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_groups_write'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.Group']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'preparation_date': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime(2012, 7, 5, 0, 0)'})
        },
        'labrack.samplecontent': {
            'Meta': {'object_name': 'SampleContent'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amount_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'mass_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentration_unit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'concentration_unit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'sample': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'sameplecontent'", 'to': "orm['labrack.Sample']"})
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