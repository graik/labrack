# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Unit'
        db.create_table('labrack_unit', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('conversion', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('unitType', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('labrack', ['Unit'])

        # Adding model 'Location'
        db.create_table('labrack_location', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='location_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('temperature', self.gf('django.db.models.fields.FloatField')()),
            ('room', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 4, 9, 0, 0))),
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
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 4, 9, 0, 0))),
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
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 4, 9, 0, 0))),
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

        # Adding model 'SelectiveMarker'
        db.create_table('labrack_selectivemarker', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['SelectiveMarker'])

        # Adding model 'Document'
        db.create_table('labrack_document', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('docfile', self.gf('django.db.models.fields.files.FileField')(max_length=100)),
        ))
        db.send_create_signal('labrack', ['Document'])

        # Adding model 'DnaSequenceAnnotation'
        db.create_table('labrack_dnasequenceannotation', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('bioStart', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('bioEnd', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('subComponent', self.gf('django.db.models.fields.related.ForeignKey')(related_name='subComponentOf', to=orm['labrack.DnaComponent'])),
            ('componentAnnotated', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='annotatedForComponent', null=True, to=orm['labrack.DnaComponent'])),
            ('strand', self.gf('django.db.models.fields.CharField')(max_length=1, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['DnaSequenceAnnotation'])

        # Adding model 'ProteinSequenceAnnotation'
        db.create_table('labrack_proteinsequenceannotation', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('bioStart', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('bioEnd', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('subComponent', self.gf('django.db.models.fields.related.ForeignKey')(related_name='subComponentOf', to=orm['labrack.ProteinComponent'])),
            ('componentAnnotated', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='annotatedForComponent', null=True, to=orm['labrack.ProteinComponent'])),
        ))
        db.send_create_signal('labrack', ['ProteinSequenceAnnotation'])

        # Adding model 'DnaComponentType'
        db.create_table('labrack_dnacomponenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('labrack', ['DnaComponentType'])

        # Adding M2M table for field subTypeOf on 'DnaComponentType'
        db.create_table('labrack_dnacomponenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_dnacomponenttype', models.ForeignKey(orm['labrack.dnacomponenttype'], null=False)),
            ('to_dnacomponenttype', models.ForeignKey(orm['labrack.dnacomponenttype'], null=False))
        ))
        db.create_unique('labrack_dnacomponenttype_subTypeOf', ['from_dnacomponenttype_id', 'to_dnacomponenttype_id'])

        # Adding model 'ProteinComponentType'
        db.create_table('labrack_proteincomponenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('labrack', ['ProteinComponentType'])

        # Adding M2M table for field subTypeOf on 'ProteinComponentType'
        db.create_table('labrack_proteincomponenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_proteincomponenttype', models.ForeignKey(orm['labrack.proteincomponenttype'], null=False)),
            ('to_proteincomponenttype', models.ForeignKey(orm['labrack.proteincomponenttype'], null=False))
        ))
        db.create_unique('labrack_proteincomponenttype_subTypeOf', ['from_proteincomponenttype_id', 'to_proteincomponenttype_id'])

        # Adding model 'ChemicalComponentType'
        db.create_table('labrack_chemicalcomponenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('labrack', ['ChemicalComponentType'])

        # Adding M2M table for field subTypeOf on 'ChemicalComponentType'
        db.create_table('labrack_chemicalcomponenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chemicalcomponenttype', models.ForeignKey(orm['labrack.chemicalcomponenttype'], null=False)),
            ('to_chemicalcomponenttype', models.ForeignKey(orm['labrack.chemicalcomponenttype'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponenttype_subTypeOf', ['from_chemicalcomponenttype_id', 'to_chemicalcomponenttype_id'])

        # Adding model 'ChassisComponentType'
        db.create_table('labrack_chassiscomponenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('labrack', ['ChassisComponentType'])

        # Adding M2M table for field subTypeOf on 'ChassisComponentType'
        db.create_table('labrack_chassiscomponenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_chassiscomponenttype', models.ForeignKey(orm['labrack.chassiscomponenttype'], null=False)),
            ('to_chassiscomponenttype', models.ForeignKey(orm['labrack.chassiscomponenttype'], null=False))
        ))
        db.create_unique('labrack_chassiscomponenttype_subTypeOf', ['from_chassiscomponenttype_id', 'to_chassiscomponenttype_id'])

        # Adding model 'PeptideComponentType'
        db.create_table('labrack_peptidecomponenttype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('labrack', ['PeptideComponentType'])

        # Adding M2M table for field subTypeOf on 'PeptideComponentType'
        db.create_table('labrack_peptidecomponenttype_subTypeOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_peptidecomponenttype', models.ForeignKey(orm['labrack.peptidecomponenttype'], null=False)),
            ('to_peptidecomponenttype', models.ForeignKey(orm['labrack.peptidecomponenttype'], null=False))
        ))
        db.create_unique('labrack_peptidecomponenttype_subTypeOf', ['from_peptidecomponenttype_id', 'to_peptidecomponenttype_id'])

        # Adding model 'Component'
        db.create_table('labrack_component', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='component_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 4, 9, 0, 0))),
            ('uri', self.gf('django.db.models.fields.URLField')(max_length=200, null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='planning', max_length=30)),
        ))
        db.send_create_signal('labrack', ['Component'])

        # Adding M2M table for field owners on 'Component'
        db.create_table('labrack_component_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_component_owners', ['component_id', 'user_id'])

        # Adding M2M table for field variantOf on 'Component'
        db.create_table('labrack_component_variantOf', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_component', models.ForeignKey(orm['labrack.component'], null=False)),
            ('to_component', models.ForeignKey(orm['labrack.component'], null=False))
        ))
        db.create_unique('labrack_component_variantOf', ['from_component_id', 'to_component_id'])

        # Adding model 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('labrack', ['ChemicalComponent'])

        # Adding M2M table for field componentType on 'ChemicalComponent'
        db.create_table('labrack_chemicalcomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chemicalcomponent', models.ForeignKey(orm['labrack.chemicalcomponent'], null=False)),
            ('chemicalcomponenttype', models.ForeignKey(orm['labrack.chemicalcomponenttype'], null=False))
        ))
        db.create_unique('labrack_chemicalcomponent_componentType', ['chemicalcomponent_id', 'chemicalcomponenttype_id'])

        # Adding model 'Chassis'
        db.create_table('labrack_chassis', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('labrack', ['Chassis'])

        # Adding M2M table for field componentType on 'Chassis'
        db.create_table('labrack_chassis_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('chassis', models.ForeignKey(orm['labrack.chassis'], null=False)),
            ('chassiscomponenttype', models.ForeignKey(orm['labrack.chassiscomponenttype'], null=False))
        ))
        db.create_unique('labrack_chassis_componentType', ['chassis_id', 'chassiscomponenttype_id'])

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

        # Adding model 'PeptideComponent'
        db.create_table('labrack_peptidecomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['PeptideComponent'])

        # Adding M2M table for field componentType on 'PeptideComponent'
        db.create_table('labrack_peptidecomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('peptidecomponent', models.ForeignKey(orm['labrack.peptidecomponent'], null=False)),
            ('peptidecomponenttype', models.ForeignKey(orm['labrack.peptidecomponenttype'], null=False))
        ))
        db.create_unique('labrack_peptidecomponent_componentType', ['peptidecomponent_id', 'peptidecomponenttype_id'])

        # Adding model 'ProteinComponent'
        db.create_table('labrack_proteincomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('GenBankfile', self.gf('django.db.models.fields.files.FileField')(max_length=100, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['ProteinComponent'])

        # Adding M2M table for field componentType on 'ProteinComponent'
        db.create_table('labrack_proteincomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('proteincomponent', models.ForeignKey(orm['labrack.proteincomponent'], null=False)),
            ('proteincomponenttype', models.ForeignKey(orm['labrack.proteincomponenttype'], null=False))
        ))
        db.create_unique('labrack_proteincomponent_componentType', ['proteincomponent_id', 'proteincomponenttype_id'])

        # Adding model 'DnaComponent'
        db.create_table('labrack_dnacomponent', (
            ('component_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Component'], unique=True, primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('translatesTo', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='encodedBy', null=True, to=orm['labrack.ProteinComponent'])),
            ('optimizedFor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['labrack.Chassis'], null=True, blank=True)),
            ('circular', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('GenBankfile', self.gf('django.db.models.fields.files.FileField')(max_length=100, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['DnaComponent'])

        # Adding M2M table for field componentType on 'DnaComponent'
        db.create_table('labrack_dnacomponent_componentType', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('dnacomponent', models.ForeignKey(orm['labrack.dnacomponent'], null=False)),
            ('dnacomponenttype', models.ForeignKey(orm['labrack.dnacomponenttype'], null=False))
        ))
        db.create_unique('labrack_dnacomponent_componentType', ['dnacomponent_id', 'dnacomponenttype_id'])

        # Adding model 'SampleCollection'
        db.create_table('labrack_samplecollection', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='samplecollection_created_by', null=True, to=orm['auth.User'])),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=25)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('labrack', ['SampleCollection'])

        # Adding M2M table for field owners on 'SampleCollection'
        db.create_table('labrack_samplecollection_owners', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('samplecollection', models.ForeignKey(orm['labrack.samplecollection'], null=False)),
            ('user', models.ForeignKey(orm['auth.user'], null=False))
        ))
        db.create_unique('labrack_samplecollection_owners', ['samplecollection_id', 'user_id'])

        # Adding model 'Sample'
        db.create_table('labrack_sample', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='sample_created_by', null=True, to=orm['auth.User'])),
            ('displayId', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('container', self.gf('django.db.models.fields.related.ForeignKey')(related_name='samples', to=orm['labrack.Container'])),
            ('aliquotNr', self.gf('django.db.models.fields.PositiveIntegerField')(null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.CharField')(default='ok', max_length=30)),
            ('reference_status', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('preparation_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 4, 9, 0, 0))),
            ('registration_date', self.gf('django.db.models.fields.DateField')(default=datetime.datetime(2013, 4, 9, 0, 0))),
            ('solvent', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
            ('concentration', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('concentrationUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='concUnit', null=True, to=orm['labrack.Unit'])),
            ('amount', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('amountUnit', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='amouUnit', null=True, to=orm['labrack.Unit'])),
            ('historyDescription', self.gf('django.db.models.fields.TextField')(max_length=255, null=True, blank=True)),
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

        # Adding M2M table for field sampleCollection on 'Sample'
        db.create_table('labrack_sample_sampleCollection', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sample', models.ForeignKey(orm['labrack.sample'], null=False)),
            ('samplecollection', models.ForeignKey(orm['labrack.samplecollection'], null=False))
        ))
        db.create_unique('labrack_sample_sampleCollection', ['sample_id', 'samplecollection_id'])

        # Adding model 'SampleProvenance'
        db.create_table('labrack_sampleprovenance', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample', self.gf('django.db.models.fields.related.ForeignKey')(related_name='sampleProvenance', to=orm['labrack.Sample'])),
            ('sample_source', self.gf('django.db.models.fields.related.ForeignKey')(related_name='sampleSource', null=True, to=orm['labrack.Sample'])),
            ('provenanceType', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['SampleProvenance'])

        # Adding model 'DnaSample'
        db.create_table('labrack_dnasample', (
            ('sample_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Sample'], unique=True, primary_key=True)),
            ('dnaConstruct', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='dnaSample', null=True, to=orm['labrack.DnaComponent'])),
            ('inChassis', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='chassisSample', null=True, to=orm['labrack.Chassis'])),
            ('derivedFrom', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='derivedFromSample', null=True, to=orm['labrack.DnaSample'])),
            ('provenanceType', self.gf('django.db.models.fields.CharField')(default='', max_length=30, null=True, blank=True)),
            ('dna_sequenced', self.gf('django.db.models.fields.CharField')(default='', max_length=30, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['DnaSample'])

        # Adding model 'ChassisSample'
        db.create_table('labrack_chassissample', (
            ('sample_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['labrack.Sample'], unique=True, primary_key=True)),
            ('chassis', self.gf('django.db.models.fields.related.ForeignKey')(related_name='Cell', to=orm['labrack.Chassis'])),
            ('derivedFrom', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='derivedFromSample', null=True, to=orm['labrack.ChassisSample'])),
            ('provenanceType', self.gf('django.db.models.fields.CharField')(default='', max_length=30, null=True, blank=True)),
        ))
        db.send_create_signal('labrack', ['ChassisSample'])


    def backwards(self, orm):
        # Removing unique constraint on 'Sample', fields ['displayId', 'container']
        db.delete_unique('labrack_sample', ['displayId', 'container_id'])

        # Deleting model 'Unit'
        db.delete_table('labrack_unit')

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

        # Deleting model 'SelectiveMarker'
        db.delete_table('labrack_selectivemarker')

        # Deleting model 'Document'
        db.delete_table('labrack_document')

        # Deleting model 'DnaSequenceAnnotation'
        db.delete_table('labrack_dnasequenceannotation')

        # Deleting model 'ProteinSequenceAnnotation'
        db.delete_table('labrack_proteinsequenceannotation')

        # Deleting model 'DnaComponentType'
        db.delete_table('labrack_dnacomponenttype')

        # Removing M2M table for field subTypeOf on 'DnaComponentType'
        db.delete_table('labrack_dnacomponenttype_subTypeOf')

        # Deleting model 'ProteinComponentType'
        db.delete_table('labrack_proteincomponenttype')

        # Removing M2M table for field subTypeOf on 'ProteinComponentType'
        db.delete_table('labrack_proteincomponenttype_subTypeOf')

        # Deleting model 'ChemicalComponentType'
        db.delete_table('labrack_chemicalcomponenttype')

        # Removing M2M table for field subTypeOf on 'ChemicalComponentType'
        db.delete_table('labrack_chemicalcomponenttype_subTypeOf')

        # Deleting model 'ChassisComponentType'
        db.delete_table('labrack_chassiscomponenttype')

        # Removing M2M table for field subTypeOf on 'ChassisComponentType'
        db.delete_table('labrack_chassiscomponenttype_subTypeOf')

        # Deleting model 'PeptideComponentType'
        db.delete_table('labrack_peptidecomponenttype')

        # Removing M2M table for field subTypeOf on 'PeptideComponentType'
        db.delete_table('labrack_peptidecomponenttype_subTypeOf')

        # Deleting model 'Component'
        db.delete_table('labrack_component')

        # Removing M2M table for field owners on 'Component'
        db.delete_table('labrack_component_owners')

        # Removing M2M table for field variantOf on 'Component'
        db.delete_table('labrack_component_variantOf')

        # Deleting model 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent')

        # Removing M2M table for field componentType on 'ChemicalComponent'
        db.delete_table('labrack_chemicalcomponent_componentType')

        # Deleting model 'Chassis'
        db.delete_table('labrack_chassis')

        # Removing M2M table for field componentType on 'Chassis'
        db.delete_table('labrack_chassis_componentType')

        # Deleting model 'Collection'
        db.delete_table('labrack_collection')

        # Removing M2M table for field components on 'Collection'
        db.delete_table('labrack_collection_components')

        # Deleting model 'PeptideComponent'
        db.delete_table('labrack_peptidecomponent')

        # Removing M2M table for field componentType on 'PeptideComponent'
        db.delete_table('labrack_peptidecomponent_componentType')

        # Deleting model 'ProteinComponent'
        db.delete_table('labrack_proteincomponent')

        # Removing M2M table for field componentType on 'ProteinComponent'
        db.delete_table('labrack_proteincomponent_componentType')

        # Deleting model 'DnaComponent'
        db.delete_table('labrack_dnacomponent')

        # Removing M2M table for field componentType on 'DnaComponent'
        db.delete_table('labrack_dnacomponent_componentType')

        # Deleting model 'SampleCollection'
        db.delete_table('labrack_samplecollection')

        # Removing M2M table for field owners on 'SampleCollection'
        db.delete_table('labrack_samplecollection_owners')

        # Deleting model 'Sample'
        db.delete_table('labrack_sample')

        # Removing M2M table for field owners on 'Sample'
        db.delete_table('labrack_sample_owners')

        # Removing M2M table for field sampleCollection on 'Sample'
        db.delete_table('labrack_sample_sampleCollection')

        # Deleting model 'SampleProvenance'
        db.delete_table('labrack_sampleprovenance')

        # Deleting model 'DnaSample'
        db.delete_table('labrack_dnasample')

        # Deleting model 'ChassisSample'
        db.delete_table('labrack_chassissample')


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
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChassisComponentType']", 'null': 'True', 'blank': 'True'}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.chassiscomponenttype': {
            'Meta': {'object_name': 'ChassisComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ChassisComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.chassissample': {
            'Meta': {'ordering': "('container', 'displayId')", 'object_name': 'ChassisSample', '_ormbases': ['labrack.Sample']},
            'chassis': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'Cell'", 'to': "orm['labrack.Chassis']"}),
            'derivedFrom': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'derivedFromSample'", 'null': 'True', 'to': "orm['labrack.ChassisSample']"}),
            'provenanceType': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '30', 'null': 'True', 'blank': 'True'}),
            'sample_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Sample']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.chemicalcomponent': {
            'Meta': {'object_name': 'ChemicalComponent', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponentType']", 'null': 'True', 'blank': 'True'}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.chemicalcomponenttype': {
            'Meta': {'object_name': 'ChemicalComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ChemicalComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
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
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'component_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'component_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 9, 0, 0)'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.Component']", 'null': 'True', 'blank': 'True'})
        },
        'labrack.container': {
            'Meta': {'object_name': 'Container'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'container_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'container_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'rack': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Rack']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 9, 0, 0)'})
        },
        'labrack.dnacomponent': {
            'GenBankfile': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'Meta': {'object_name': 'DnaComponent', '_ormbases': ['labrack.Component']},
            'circular': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.DnaComponentType']", 'null': 'True', 'blank': 'True'}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['labrack.Chassis']", 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"})
        },
        'labrack.dnacomponenttype': {
            'Meta': {'object_name': 'DnaComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
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
            'sample_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Sample']", 'unique': 'True', 'primary_key': 'True'})
        },
        'labrack.dnasequenceannotation': {
            'Meta': {'object_name': 'DnaSequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'componentAnnotated': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'annotatedForComponent'", 'null': 'True', 'to': "orm['labrack.DnaComponent']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1', 'null': 'True', 'blank': 'True'}),
            'subComponent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'subComponentOf'", 'to': "orm['labrack.DnaComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.document': {
            'Meta': {'object_name': 'Document'},
            'docfile': ('django.db.models.fields.files.FileField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        'labrack.location': {
            'Meta': {'object_name': 'Location'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'location_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'location_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 9, 0, 0)'}),
            'room': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'temperature': ('django.db.models.fields.FloatField', [], {})
        },
        'labrack.peptidecomponent': {
            'Meta': {'object_name': 'PeptideComponent', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.PeptideComponentType']", 'null': 'True', 'blank': 'True'}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'})
        },
        'labrack.peptidecomponenttype': {
            'Meta': {'object_name': 'PeptideComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.PeptideComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponent': {
            'GenBankfile': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'Meta': {'object_name': 'ProteinComponent', '_ormbases': ['labrack.Component']},
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.ProteinComponentType']", 'null': 'True', 'blank': 'True'}),
            'component_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['labrack.Component']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'})
        },
        'labrack.proteincomponenttype': {
            'Meta': {'object_name': 'ProteinComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['labrack.ProteinComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'labrack.proteinsequenceannotation': {
            'Meta': {'object_name': 'ProteinSequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'componentAnnotated': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'annotatedForComponent'", 'null': 'True', 'to': "orm['labrack.ProteinComponent']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'subComponent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'subComponentOf'", 'to': "orm['labrack.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
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
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 9, 0, 0)'})
        },
        'labrack.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'aliquotNr': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'amountUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'amouUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'concentrationUnit': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'concUnit'", 'null': 'True', 'to': "orm['labrack.Unit']"}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['labrack.Container']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sample_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'historyDescription': ('django.db.models.fields.TextField', [], {'max_length': '255', 'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'sample_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"}),
            'preparation_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 9, 0, 0)'}),
            'reference_status': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'registration_date': ('django.db.models.fields.DateField', [], {'default': 'datetime.datetime(2013, 4, 9, 0, 0)'}),
            'sampleCollection': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['labrack.SampleCollection']", 'null': 'True', 'blank': 'True'}),
            'solvent': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'ok'", 'max_length': '30'})
        },
        'labrack.samplecollection': {
            'Meta': {'object_name': 'SampleCollection'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'samplecollection_created_by'", 'null': 'True', 'to': "orm['auth.User']"}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '25'}),
            'owners': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'samplecollection_owners'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['auth.User']"})
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
        'labrack.unit': {
            'Meta': {'object_name': 'Unit'},
            'conversion': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'unitType': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        }
    }

    complete_apps = ['labrack']