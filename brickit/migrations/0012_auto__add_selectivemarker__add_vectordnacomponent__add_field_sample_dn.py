# encoding: utf-8
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models

class Migration(SchemaMigration):

    def forwards(self, orm):
        
        # Adding model 'SelectiveMarker'
        db.create_table('brickit_selectivemarker', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
            ('shortDescription', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('description', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal('brickit', ['SelectiveMarker'])

        # Adding model 'VectorDnaComponent'
        db.create_table('brickit_vectordnacomponent', (
            ('dnacomponent_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['brickit.DnaComponent'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('brickit', ['VectorDnaComponent'])

        # Adding M2M table for field marker on 'VectorDnaComponent'
        db.create_table('brickit_vectordnacomponent_marker', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('vectordnacomponent', models.ForeignKey(orm['brickit.vectordnacomponent'], null=False)),
            ('selectivemarker', models.ForeignKey(orm['brickit.selectivemarker'], null=False))
        ))
        db.create_unique('brickit_vectordnacomponent_marker', ['vectordnacomponent_id', 'selectivemarker_id'])

        # Adding field 'Sample.dna'
        db.add_column('brickit_sample', 'dna', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='dna_samples', null=True, to=orm['brickit.DnaComponent']), keep_default=False)

        # Adding field 'Sample.vector'
        db.add_column('brickit_sample', 'vector', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='vector_samples', null=True, to=orm['brickit.VectorDnaComponent']), keep_default=False)

        # Adding field 'Sample.cell'
        db.add_column('brickit_sample', 'cell', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='chassis_samples', null=True, to=orm['brickit.Chassis']), keep_default=False)

        # Adding field 'Sample.protein'
        db.add_column('brickit_sample', 'protein', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='protein_samples', null=True, to=orm['brickit.ProteinComponent']), keep_default=False)


    def backwards(self, orm):
        
        # Deleting model 'SelectiveMarker'
        db.delete_table('brickit_selectivemarker')

        # Deleting model 'VectorDnaComponent'
        db.delete_table('brickit_vectordnacomponent')

        # Removing M2M table for field marker on 'VectorDnaComponent'
        db.delete_table('brickit_vectordnacomponent_marker')

        # Deleting field 'Sample.dna'
        db.delete_column('brickit_sample', 'dna_id')

        # Deleting field 'Sample.vector'
        db.delete_column('brickit_sample', 'vector_id')

        # Deleting field 'Sample.cell'
        db.delete_column('brickit_sample', 'cell_id')

        # Deleting field 'Sample.protein'
        db.delete_column('brickit_sample', 'protein_id')


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
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.Chassis']", 'null': 'True', 'blank': 'True'})
        },
        'brickit.collection': {
            'Meta': {'object_name': 'Collection'},
            'chassis': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.Chassis']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'dnaComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.DnaComponent']", 'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'proteinComponents': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ProteinComponent']", 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'brickit.componenttype': {
            'Meta': {'object_name': 'ComponentType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'subTypeOf': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'subTypes'", 'null': 'True', 'symmetrical': 'False', 'to': "orm['brickit.ComponentType']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200'})
        },
        'brickit.dnacomponent': {
            'Meta': {'object_name': 'DnaComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'optimizedFor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['brickit.Chassis']", 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'translatesTo': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'encodedBy'", 'null': 'True', 'to': "orm['brickit.ProteinComponent']"}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.DnaComponent']", 'null': 'True', 'blank': 'True'})
        },
        'brickit.proteincomponent': {
            'Meta': {'object_name': 'ProteinComponent'},
            'abstract': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'annotations': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'componentType': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ComponentType']", 'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'planning'", 'max_length': '30'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'variantOf': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.ProteinComponent']", 'null': 'True', 'blank': 'True'})
        },
        'brickit.sample': {
            'Meta': {'ordering': "('container', 'displayId')", 'unique_together': "(('displayId', 'container'),)", 'object_name': 'Sample'},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'chassis_samples'", 'null': 'True', 'to': "orm['brickit.Chassis']"}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '50', 'max_digits': '6', 'decimal_places': '2'}),
            'concentrationUnit': ('django.db.models.fields.CharField', [], {'default': "'mg/l'", 'max_length': '10'}),
            'container': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'samples'", 'to': "orm['brickit.SampleContainer']"}),
            'created': ('django.db.models.fields.DateField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'dna': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'dna_samples'", 'null': 'True', 'to': "orm['brickit.DnaComponent']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'protein_samples'", 'null': 'True', 'to': "orm['brickit.ProteinComponent']"}),
            'sampleType': ('django.db.models.fields.CharField', [], {'default': "'dna'", 'max_length': '100'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.User']", 'db_index': 'True', 'symmetrical': 'False'}),
            'vector': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'vector_samples'", 'null': 'True', 'to': "orm['brickit.VectorDnaComponent']"}),
            'vesselType': ('django.db.models.fields.CharField', [], {'default': 'None', 'max_length': '30', 'null': 'True', 'blank': 'True'})
        },
        'brickit.samplecontainer': {
            'Meta': {'ordering': "('displayId',)", 'object_name': 'SampleContainer'},
            'containerType': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'displayId': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['auth.User']", 'null': 'True', 'db_index': 'True', 'blank': 'True'})
        },
        'brickit.selectivemarker': {
            'Meta': {'ordering': "('name',)", 'object_name': 'SelectiveMarker'},
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'shortDescription': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        },
        'brickit.sequenceannotation': {
            'Meta': {'object_name': 'SequenceAnnotation'},
            'bioEnd': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'bioStart': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'precedes': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.SequenceAnnotation']", 'null': 'True', 'blank': 'True'}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1', 'null': 'True', 'blank': 'True'}),
            'uri': ('django.db.models.fields.URLField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
        },
        'brickit.vectordnacomponent': {
            'Meta': {'object_name': 'VectorDnaComponent', '_ormbases': ['brickit.DnaComponent']},
            'dnacomponent_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['brickit.DnaComponent']", 'unique': 'True', 'primary_key': 'True'}),
            'marker': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': "orm['brickit.SelectiveMarker']", 'null': 'True', 'blank': 'True'})
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
