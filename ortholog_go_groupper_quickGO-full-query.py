
# EGGNOG GROUPPER v1.5
#
# What this file do?
# This file converts downloaded id_mapping files (from data/uniprot folder) to a reduced file size. Filtering out STRING and eggNOG db related lines.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import csv
import sys
import operator

csv.field_size_limit(sys.maxsize)
print_details = True

# List of GO_SLIM terms
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
    'GO:0000278': 'mitotic cell cycle'
    }

list_of_selected_10 = {
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

# Taxon specific collectors
taxon_list = ['3702', '6239', '7227', '7955', '9606', '559292', '284812']
taxon_dict = {'3702': [], '6239': [], '7227': [], '7955': [], '9606': [], '559292': [], '284812': []}
taxon_dict_names = {
    '284812': 'S. pombe',
    '3702': 'A. thaliana',
    '4896': 'S. pombe',
    '559292': 'S. cerevisiae',
    '6239': 'C. elegans',
    '63221': 'H. sapiens', # Homo sapiens neanderthalensis (subspecies).
    '7227': 'D. melanogaster',
    '741158': 'H. sapiens', # Homo sapiens ssp. Denisova (subspecies).
    '7955': 'D. rerio',
    '9606': 'H. sapiens'
}

###############################
# Reading QuickGO export file #
###############################

# file for GO_SLIM generic
filename_go = "data/go/QuickGO-annotations-GOslim-generic-20220126.tsv"
# file for our selected 10
# filename_go = "data/go/QuickGO-annotations-special-10-20220127.tsv"
# list_of_goSLIM_generic = list_of_selected_10

go_goslim = {}
for go_id in list_of_goSLIM_generic.keys():
    go_goslim[go_id] = []

go_subgo = {}
with open(filename_go, newline = '') as f:
    reader = csv.DictReader(f, fieldnames = ('type', 'uniprot', 'name', 'qualifier', 'goslim', 'subgo', 'longname', 'eco_id', 'go_evidence', 'reference', 'with_from', 'taxid', 'assigned_by', 'annot_ext', 'go_aspect'), delimiter = '\t')
    counter = 0
    counter_for_goslim = 0
    try:
        for row in reader:
            if row['type'] != 'UniProtKB':
                continue
            if row['goslim'] not in list_of_goSLIM_generic.keys():
                print(row['subgo'], row['goslim'])
                counter_for_goslim += 1
            go_goslim[row['goslim']].append(row['uniprot'])
            go_subgo[row['uniprot']] = row['subgo']
            counter += 1
            # if counter == 100:
            #     break
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

print(f"Parser have {counter} lines processed from {filename_go} file.")

# In case our selection of GOSlim contains more origins of GOs - print out the number
if counter_for_goslim > 0:
    print(counter_for_goslim)

################################
# Reading UniProt Convert file #
################################

# Source filename for eggNOG_data (based on Uniprot database)
filename_uniprot_convert = "data/uniprot_convert_merged.tsv"
# Declare eggNOG database output variable
eggNOG_database = {}

with open(filename_uniprot_convert, newline = '') as f:
    reader = csv.DictReader(f, fieldnames = ('uniprot', 'db', 'taxid', 'groupid'), delimiter = '\t')
    counter = 0
    try:
        for row in reader:
            counter += 1

            # Filtering eggNOG DB out
            if row['db'] == 'eggNOG':

                if row['groupid'] in eggNOG_database:
                    if row['taxid'] in eggNOG_database[row['groupid']]:
                        eggNOG_database[row['groupid']][row['taxid']].append(row['uniprot'])
                    else:
                        eggNOG_database[row['groupid']][row['taxid']] = [row['uniprot']]
                else:
                    eggNOG_database[row['groupid']] = {row['taxid']: [row['uniprot']]}

    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

    print("Parser have", counter, "lines processed from eggNOG merged file.")


################################
# Groupping #
# ################################

counter_for_slim = 0
for this_goslim in list_of_goSLIM_generic.keys():

    counter_for_slim += 1

    write_the_output = []
    write_the_output_hit = []
    write_the_output_hit_dict_4items = {}
    write_the_output_hit_dict_4order = {}

    counter = 0
    for groupid in eggNOG_database:

        this_line = ""
        this_mezok = {}
        this_hit_spec_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}
        this_hit_prot_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}
        this_total_spec_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}
        this_total_prot_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}

        for sor_taxid in taxon_list:
            sor_taxid2 = "H" + sor_taxid + "H"
            if sor_taxid in eggNOG_database[groupid]:
                this_mezok[sor_taxid] = []
                this_mezok[sor_taxid2] = []
                this_total_spec_count[sor_taxid] = 1
                for uniprot in eggNOG_database[groupid][sor_taxid]:

                    this_mezok[sor_taxid].append(uniprot)
                    this_total_prot_count[sor_taxid] += 1
                    if uniprot in go_goslim[this_goslim]:
                        this_hit_spec_count[sor_taxid] = 1
                        this_hit_prot_count[sor_taxid] += 1
                        this_mezok[sor_taxid2].append(uniprot)
            else:
                this_mezok[sor_taxid] = []
                this_mezok[sor_taxid2] = []

        average_hm_list = []

        hit_spec_count = 0
        hit_prot_count = 0
        total_spec_count = 0
        total_prot_count = 0
        for taxid in taxon_list:
            hit_spec_count += this_hit_spec_count[taxid]
            hit_prot_count += this_hit_prot_count[taxid]
            total_spec_count += this_total_spec_count[taxid]
            total_prot_count += this_total_prot_count[taxid]

            if this_hit_spec_count[taxid] > 0:
                average_hm_list.append(float(this_hit_prot_count[taxid] / this_total_prot_count[taxid]))

        average_hm_total = 0
        if len(average_hm_list) == 0:
            average_hm = 0
        else:
            for value in average_hm_list:
                average_hm_total += value
            average_hm = float(average_hm_total / total_spec_count)

        if hit_spec_count > 0:
            total_hm = hit_prot_count / total_prot_count
        else:
            total_hm = 0

        # Writing out the cells in the order of appearance

        this_line += groupid + "\t"
        this_line += str(float("{:.5f}".format(average_hm))) + "\t"
        this_line += str(float("{:.5f}".format(total_hm))) + "\t"
        this_line += str(hit_prot_count) + "\t"
        this_line += str(total_prot_count) + "\t"
        this_line += str(hit_spec_count) + "\t"
        this_line += str(total_spec_count) + "\t"

        if print_details:
            # Adding detailed information
            for taxid in taxon_list:
                this_line += str(this_hit_spec_count[taxid]) + "\t"
                this_line += str(this_hit_prot_count[taxid]) + "\t"
                #this_line += str(this_total_spec_count[taxid]) + "\t"
                #this_line += str(this_total_prot_count[taxid]) + "\t"

                if this_total_prot_count[taxid] > 0:
                    fraction = this_hit_prot_count[taxid] / this_total_prot_count[taxid]
                else:
                    fraction = 0

                #this_line += str(float("{:.5f}".format(fraction))) + "\t"
        else:

            for sor_taxid in taxon_list:
                this_line += ",".join(this_mezok[sor_taxid]) + "\t"
            for sor_taxid in taxon_list:
                sor_taxid2 = "H" + sor_taxid + "H"
                this_line += ",".join(this_mezok[sor_taxid2]) + "\t"

        counter += 1
        #write_the_output.append(this_line)
        if hit_spec_count > 0:
            if total_spec_count > 3:
                write_the_output_hit.append(this_line)
                write_the_output_hit_dict_4items[counter] = this_line
                write_the_output_hit_dict_4order[counter] = int(average_hm * 10000000)

    print(f"Cycle Nr {counter_for_slim} - Parser have {counter} elements processed from {this_goslim}.")

    sorted_dicdata = sorted(write_the_output_hit_dict_4order.items(), key = operator.itemgetter(1), reverse = True)

    if print_details:
        export_filename = "data/go-" + this_goslim.replace('GO:', '') + "-ordered-detailed-2022-01-correct.tsv"
    else:
        export_filename = "data/go-" + this_goslim.replace('GO:', '') + "-ordered-2022-01.tsv"

    counter = 0

    with open(export_filename, mode = 'w') as export_file:
        writer = csv.writer(export_file, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)

        this_line = "Group ID\t"
        this_line += "Average H/M\t"
        this_line += "Total H/M\t"
        this_line += "Hit Proteins\t"
        this_line += "Total Proteins\t"
        this_line += "Hit Species\t"
        this_line += "Total Species\t"

        for sor_taxid in taxon_list:
            this_line += taxon_dict_names[sor_taxid] + " Total\t"
        for sor_taxid in taxon_list:
            this_line += taxon_dict_names[sor_taxid] + " Hit\t"

        this_line = this_line.split("\t")
        writer.writerow(this_line)

        for rowid in sorted_dicdata:
            counter += 1
            this_line = write_the_output_hit_dict_4items[rowid[0]].split("\t")
            writer.writerow(this_line)

        print(f"Cycle Nr {counter_for_slim} - Parser have {counter} lines wrote in {export_filename} (merged) file.")
