# Ortholog Parser / STRING DB Reader
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins in a given bottle.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#
# File for functions

# Import libraries
import csv
import sys
import random
import logging
import pandas as pd

file_dir = "data/string-values/202202/"
file_source = "pvalues_go-0000902_detailed_20210617-062529-863990.tsv"
how_together = 6
export_excel_filename = str(how_together) + "together_202202.xlsx"

# sources = {"GO:0000502": "pvalues_go-0000502_detailed_20210527-140415-591244.tsv",
#            "GO:0000902": "pvalues_go-0000902_detailed_20210617-062529-863990.tsv",
#            "GO:0000910": "pvalues_go-0000910_detailed_20210617-063520-677966.tsv",
#            "GO:0002376": "pvalues_go-0002376_detailed_20210617-065046-555579.tsv",
#            "GO:0003013": "pvalues_go-0003013_detailed_20210617-071021-229058.tsv",
#            "GO:0005975": "pvalues_go-0005975_detailed_20210617-071333-249858.tsv",
#            "GO:0006099": "pvalues_go-0006099_detailed_20210527-170952-549387.tsv",
#            "GO:0006259": "pvalues_go-0006259_detailed_20210706-144222-468228.tsv",
#            "GO:0006397": "pvalues_go-0006397_detailed_20210706-150214-942212.tsv",
#            "GO:0006399": "pvalues_go-0006399_detailed_20210706-153436-644640.tsv",
#            "GO:0006412": "pvalues_go-0006412_detailed_20210520-164627-847050.tsv",
#            "GO:0006629": "pvalues_go-0006629_detailed_20210520-160012-631656.tsv",
#            "GO:0006914": "pvalues_go-0006914_detailed_20210520-163643-811997.tsv",
#            "GO:0007049": "pvalues_go-0007049_detailed_20210520-165833-216219.tsv",
#            "GO:0007568": "pvalues_go-0007568_detailed_20210520-171115-983076.tsv",
#            "GO:0008361": "pvalues_go-0008361_detailed_20210520-172247-251322.tsv",
#            "GO:0009295": "pvalues_go-0009295_detailed_20210527-153107-149308.tsv",
#            "GO:0051301": "pvalues_go-0051301_detailed_20210520-172831-430456.tsv",
#            "GO:0051726": "pvalues_go-0051726_detailed_20210504-154234-473575.tsv"}

list_of_goSLIM_generic = {
	'GO:1901135': 'carbohydrate derivative metabolic process',
	'GO:0140053': 'mitochondrial gene expression',
	'GO:0140014': 'mitotic nuclear division',
	'GO:0140013': 'meiotic nuclear division',
	'GO:0098754': 'detoxification',
	'GO:0098542': 'defense response to other organism',
	'GO:0072659': 'protein localization to plasma membrane',
	'GO:0071941': 'nitrogen cycle metabolic process',
	'GO:0071554': 'cell wall organization or biogenesis',
	'GO:0065003': 'protein-containing complex assembly',
	'GO:0061024': 'membrane organization',
	# 'GO:0061007': 'hepaticobiliary system process',
	'GO:0055086': 'nucleobase-containing small molecule metabolic process',
	'GO:0055085': 'transmembrane transport',
	'GO:0055065': 'metal ion homeostasis',
	'GO:0051604': 'protein maturation',
	'GO:0050886': 'endocrine process',
	'GO:0050877': 'nervous system process',
	'GO:0048870': 'cell motility',
	'GO:0048856': 'anatomical structure development',
	'GO:0044782': 'cilium organization',
	'GO:0042254': 'ribosome biogenesis',
	'GO:0042060': 'wound healing',
	'GO:0036211': 'protein modification process',
	'GO:0034330': 'cell junction organization',
	'GO:0032200': 'telomere organization',
	'GO:0031047': 'gene silencing by RNA',
	'GO:0030198': 'extracellular matrix organization',
	'GO:0030163': 'protein catabolic process',
	'GO:0030154': 'cell differentiation',
	'GO:0023052': 'signaling',
	'GO:0022600': 'digestive system process',
	'GO:0022414': 'reproductive process',
	'GO:0016192': 'vesicle-mediated transport',
	'GO:0016073': 'snRNA metabolic process',
	'GO:0016071': 'mRNA metabolic process',
	'GO:0015979': 'photosynthesis',
	'GO:0012501': 'programmed cell death',
	'GO:0007568': 'aging',
	'GO:0007163': 'establishment or maintenance of cell polarity',
	'GO:0007155': 'cell adhesion',
	'GO:0007059': 'chromosome segregation',
	'GO:0007040': 'lysosome organization',
	'GO:0007031': 'peroxisome organization',
	'GO:0007018': 'microtubule-based movement',
	'GO:0007010': 'cytoskeleton organization',
	'GO:0007005': 'mitochondrion organization',
	'GO:0006954': 'inflammatory response',
	'GO:0006914': 'autophagy',
	'GO:0006913': 'nucleocytoplasmic transport',
	'GO:0006886': 'intracellular protein transport',
	'GO:0006790': 'sulfur compound metabolic process',
	'GO:0006766': 'vitamin metabolic process',
	'GO:0006629': 'lipid metabolic process',
	'GO:0006575': 'cellular modified amino acid metabolic process',
	'GO:0006520': 'cellular amino acid metabolic process',
	'GO:0006486': 'protein glycosylation',
	'GO:0006457': 'protein folding',
	'GO:0006399': 'tRNA metabolic process',
	'GO:0006355': 'regulation of transcription, DNA-templated',
	'GO:0006351': 'transcription, DNA-templated',
	'GO:0006325': 'chromatin organization',
	'GO:0006310': 'DNA recombination',
	'GO:0006281': 'DNA repair',
	'GO:0006260': 'DNA replication',
	'GO:0006091': 'generation of precursor metabolites and energy',
	'GO:0005975': 'carbohydrate metabolic process',
	'GO:0003016': 'respiratory system process',
	'GO:0003014': 'renal system process',
	'GO:0003013': 'circulatory system process',
	'GO:0003012': 'muscle system process',
	'GO:0002376': 'immune system process',
	'GO:0002181': 'cytoplasmic translation',
	'GO:0000910': 'cytokinesis',
	'GO:0000278': 'mitotic cell cycle',

	# NEW LINER: selected 10 proteiens

	'GO:0007049': 'cell cycle',
	'GO:0000902': 'cell morphogenesis',
	'GO:0006259': 'DNA metabolic process',
	'GO:0008361': 'regulation of cell size',
	'GO:0051726': 'regulation of cell cycle',
	'GO:0051301': 'cell division',
	'GO:0006412': 'translation',
	'GO:0006099': 'tricarboxylic acid cycle',
	'GO:0000502': 'proteasome complex',
	'GO:0009295': 'nucleoid'
}

go_ids_list = []
for go_id in list_of_goSLIM_generic.keys():
	go_ids_list.append(go_id.replace("GO:", "GO_"))

# go_ids_list = ["GO_0008361", "GO_0009295", "GO_0002376", "GO_0006399", "GO_0000902", "GO_0000910",
#            "GO_0003013", "GO_0006099", "GO_0005975", "GO_0006259", "GO_0006914", "GO_0006629",
#            "GO_0051726", "GO_0051301", "GO_0006397", "GO_0007568", "GO_0007049", "GO_0000502",
#            "GO_0006412"]

sources = {}

for goid in go_ids_list:
	sources[goid.replace("_", ":")] = "qucikgo_export_" + goid.replace("GO_", "GO-") + "_" + str(how_together) + "together_202202.tsv"

# sources = {"GO:0009295": "pvalues_go-0009295_detailed_20210527-153107-149308.tsv"}

bottles_list1 = (0.00, 0.33, 0.66)
bottles_list2 = (0.33, 0.66, 1.00)
taxon_list = (9606, 7955, 6239, 3702, 7227, 4896, 4932)
# corr_types_list = ('pearson', 'kendall', 'spearman')
corr_types_list = ('pearson', 'spearman')

export_list = []
df_dict = {}

for goid, filename in sources.items():

	df = pd.read_csv(file_dir + filename,
	                 sep = '\t',
	                 names = ["taxid",
	                          "bottle_1",
	                          "bottle_2",
	                          "hm_average",
	                          "number_of_nodes",
	                          "number_of_edges",
	                          "average_node_degree",
	                          "local_clustering_coefficient",
	                          "expected_number_of_edges",
	                          "p_value"]
	                 )

	for corr_type in corr_types_list:

		for taxon in taxon_list:

			corr = df[df["taxid"] == taxon].corr(method = corr_type)
			export_row = [goid, corr_type]
			for element in corr["hm_average"].iteritems():
				if element[0] == "taxid":
					export_row.append(taxon)
				elif element[0] == "bottle_1":
					export_row.append("ALL")
				elif element[0] == "bottle_2":
					export_row.append("ALL")
				else:
					export_row.append(element[1])
			export_list.append(export_row)

			bottle_num = 0

			# for bottle_1 in bottles_list1:
			#
			# 	corr = df[(df["taxid"] == taxon) & (df["bottle_1"] == bottle_1)].corr(method=corr_type)
			# 	export_row = [goid, corr_type]
			# 	for element in corr["hm_average"].iteritems():
			# 		if element[0] == "taxid":
			# 			export_row.append(taxon)
			# 		elif element[0] == "bottle_1":
			# 			export_row.append(bottle_1)
			# 		elif element[0] == "bottle_2":
			# 			export_row.append(bottles_list2[bottle_num])
			# 		else:
			# 			export_row.append(element[1])
			# 	export_list.append(export_row)
			#
			# 	bottle_num += 1

	# Rprint(export_list)

	dfx = pd.DataFrame(export_list)
	dfx.columns = ["goid",
	               "corr_type",
	               "taxid",
	               "bottle_1",
	               "bottle_2",
	               "hm_average",
	               "number_of_nodes",
	               "number_of_edges",
	               "average_node_degree",
	               "local_clustering_coefficient",
	               "expected_number_of_edges",
	               "p_value"]

	df_dict[goid] = dfx
	export_list = []

with pd.ExcelWriter(export_excel_filename) as writer:
	for goid, df_item in df_dict.items():
		df_item.to_excel(writer, sheet_name = goid.replace(":", "_"))
