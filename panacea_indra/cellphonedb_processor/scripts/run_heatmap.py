import os
import sys

wd = '/Users/sbunga/gitHub/panacea_indra/cellphone_db/output'

samples = {
		'zymo': (
			os.path.join(wd, 'cellphonedb_output/neuro_zymo/pvalues.txt'), 
			os.path.join(wd, 'meta_data/neuro_Zymo_meta_table.tsv')
			),
		'incision': (
			os.path.join(wd, 'cellphonedb_output/neuro_incision/pvalues.txt'), 
			os.path.join(wd, 'meta_data/neuro_Incision_meta_table.tsv')
			),
		'healthy': (
			os.path.join(wd, 'cellphonedb_output/neuro_healthy/pvalues.txt'), 
			os.path.join(wd, 'meta_data/neuro_Healthy_meta_table.tsv')
			),
		'saline': (
			os.path.join(wd, 'cellphonedb_output/neuro_saline/pvalues.txt'), 
			os.path.join(wd, 'meta_data/neuro_Saline_meta_table.tsv')
			),
		'uvb': (
			os.path.join(wd, 'cellphonedb_output/neuro_uvb/pvalues.txt'), 
			os.path.join(wd, 'meta_data/neuro_UVB_meta_table.tsv')
			),
		'sham': (
			os.path.join(wd, 'cellphonedb_output/neuro_sham/pvalues.txt'), 
			os.path.join(wd, 'meta_data/neuro_Sham_meta_table.tsv')
			)


}


for k in samples.keys():
      
    os.system('cellphonedb plot heatmap_plot '+samples[k][1] + ' --pvalues-path '+ samples[k][0] + \
    	' --output-path '+os.path.join(os.path.dirname(samples[k][0]), 'plots'))







