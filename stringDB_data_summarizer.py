
# EGGNOG RANDOM SELECTOR
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins, randomly selecting the pool.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

# Import libraries
import pandas as pd
# Import local functions
from stringDB_functions import *
# Import local variables
from stringDB_variables import *

#######################
# SET FILE PARAMETERS #
#######################

isTest = False
num_cycles = 40
num_request_per_cycle = 10
num_proteins = 6
dir_export = "export/"
dir_log = "logs/"
filename_ext = "-2022-01"

# list_goids = ["go-0008361", "go-0002376", "go-0009295", "go-0000902",
#               "go-0006099", "go-0003013", "go-0000502", "go-0006399",
#               "go-0000910", "go-0005975"]

# list_goids = ["go-0006629", "go-0006914", "go-0007568", "go-0006259",
#               "go-0051726", "go-0051301", "go-0006397", "go-0007049",
#               "go-0006412"]

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
    'GO:0061007': 'hepaticobiliary system process',
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

# list_of_goSLIM_generic = list_of_selected_10

list_goids = []
for go_id in list_of_goSLIM_generic.keys():
	list_goids.append(go_id.replace("GO:", "go-"))

for str_goid in list_goids:

	#str_goid = "go-0008361"
	log_filename1 = dir_log + "data_" + str_goid + "_" + str(num_proteins) + "together_general_" + current_time_abbrev + ".tsv"
	log_filename2 = dir_log + "data_" + str_goid + "_" + str(num_proteins) + "together_detailed_" + current_time_abbrev + ".tsv"

	# START LOGGING
	logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename1, level=logging.DEBUG)

	# Counter of Requests
	counter_requests = 0

	########################################################
	# Open GO file in Pandas as a DataFrame for this TaxID #
	########################################################
	selected_cols = [	"Group ID",
	                     "Average H/M",
	                     "Total H/M",
	                     "Hit Proteins",
	                     "Total Proteins",
	                     "Hit Species",
	                     "Total Species",
	                     "A. thaliana Total",
	                     "C. elegans Total",
	                     "D. melanogaster Total",
	                     "D. rerio Total",
	                     "H. sapiens Total",
	                     "S. cerevisiae Total",
	                     "S. pombe Total" ]

	column_name_hm_value = "Average H/M"

	go = pd.read_csv("data/" + str_goid + "-ordered" + filename_ext + ".tsv", sep="\t", usecols=selected_cols)

	# Print & Log pd read
	print(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")
	logging.info(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")

	# Logging
	log_calls = []

	# Create P-value array for the whole taxid
	p_val_array = {}

	num_all_groups = len(go)-1


	break
