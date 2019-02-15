# prepare_model.py

# the function make_string_labels_independent may cause a Open Babel error. This is expected and should not cause adverse effects

import cantera_tools as ctt
from rmgpy.chemkin import loadChemkinFile
from rmgpy.tools.isotopes import cluster
import os
import analysis_methods as am

## mainPathcontains the file with name `chem.inp` and `species_dictionary.txt`. 
mainPath = '../mechanisms/drg_model'

# make cti mechanism from inp file
#ctt.obtain_cti_file_nicely_named(mainPath,original_ck_file='chem.inp')

# create the file mapping isotopologues to isotopes

chemkinPath = os.path.join(mainPath, 'chem.inp')
speciesDictPath = os.path.join(mainPath,'species_dictionary.txt')
species, reactions = loadChemkinFile(chemkinPath, speciesDictPath, readComments = False, useChemkinNames=False)
ctt.make_string_labels_independent(species)

spec_clusters = cluster(species)

cluster_info = am.getIsotomoperInfo(spec_clusters)
cluster_info['not_r_enriched'] = cluster_info['enriched_atoms'] - cluster_info['r_enriched']
cluster_info['not_r_unenriched'] = cluster_info['unenriched_atoms'] - cluster_info['r_unenriched']
cluster_info.to_csv(os.path.join(mainPath, 'isotopomer_cluster_info.csv'),index_label='name')
