
# EGGNOG RANDOM SELECTOR
#
# What this file do?
# x
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import os
import csv
import sys
import random
import requests ## python -m pip install requests
import pandas as pd
from datetime import datetime

csv.field_size_limit(sys.maxsize)

#
# SET FILE
#

num_cycles = 10
num_request_per_cycle = 10
num_bottles = 5
str_goid = "go-0051726"


# Creating timestamp for output filename
now = datetime.now()
current_time_abbrev = now.strftime("%Y%m%d-%H%M%S-%f")

#
# Taxon files
#

taxon_list = ('9606', '7955', '6239', '3702', '7227', '4896', '4932')
taxon_dict = {'9606': 'H. sapiens', '7955': 'D. rerio', '6239': 'C. elegans', '3702': 'A. thaliana', '7227': 'D. melanogaster', '4896': 'S. pombe', '4932': 'S. cerevisiae'}
taxon_dict_go = {'9606': 'H. sapiens Hit', '7955': 'D. rerio Hit', '6239': 'C. elegans Hit', '3702': 'A. thaliana Hit', '7227': 'D. melanogaster Hit', '4896': 'S. pombe Hit', '4932': 'S. cerevisiae Hit'}

taxid = taxon_list[0]

#
# Read UniProt file
#

filename = "data/uniprot_convert_" + taxid + ".tsv"
write_lines_all = {}
convert_uniprot = {}
name_uniprot = {}
list_of_uniprots = []

with open(filename, newline='') as f:
	reader = csv.DictReader(f, fieldnames=('uniprot', 'db', 'taxid', 'selector'), delimiter='\t')
	counter = 0
	try:
		for row in reader:

			# Filtering STRING to convert Uniprot to STRING db id
			if row['db'] == 'convert':
				if row['uniprot'] in convert_uniprot:
					continue
				else:
					convert_uniprot[row['uniprot']] = row['selector']
					list_of_uniprots.append(row['uniprot'])
					counter += 1

			# Filtering STRING convert out
			if row['db'] == 'Gene_Name':
				if row['uniprot'] in name_uniprot:
					continue
				else:
					name_uniprot[row['uniprot']] = row['selector']

	except csv.Error as e:
		sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

	print(f"FileReader have {counter} lines load from {taxid} list/file.")

#
# Design Bottles
#

bottles = {1.0: [], 0.8: [], 0.6: [], 0.4: [], 0.2: []}
if num_bottles == 5:
	bottles_list = (0, 0.2, 0.4, 0.6, 0.8, 1)
elif num_bottles == 10:
	bottles_list = (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

#
# Open GO file in Pandas as a DataFrame
#

selected_cols = ["Group ID", "Average H/M", "Total H/M", "Hit Proteins", "Total Proteins", "Hit Species", "Total Species"]
this_taxon_column_name = taxon_dict_go[taxon_list[0]]
selected_cols.append(this_taxon_column_name)

go = pd.read_csv("data/" + str_goid + "-ordered.tsv", sep="\t", usecols=selected_cols)

# Logging
log_calls = []

bottle_counter = 0
for num in bottles_list:

	if bottle_counter == (len(bottles_list)-1):
		continue

	nr1 = bottles_list[bottle_counter]
	nr2 = bottles_list[bottle_counter + 1]
	go2 = go["Average H/M"].between(nr1, nr2)

	bottle_counter += 1

	print(f"Bottle between {nr1}-{nr2} started.")

	# Retrieving a protein list from the filtered pandas dataset of this species
	list_of_bottle_proteins = []
	for prot in go[go2][this_taxon_column_name]:
		prot = str(prot)

		# When there is no hit protein for this species: SKIP
		if prot == "nan":
			continue
		# When there are more than one hit protein, randomly selection one
		if (prot.find(',') != -1):
			this_prots = prot.split(",")
			selected = random.choice(this_prots)
		else:
			selected = prot

		# Adding proteins to a list
		list_of_bottle_proteins.append(selected)

	for k in range(1, (num_cycles+1)):
		responses_array = []
		p_values = []

		print(f"STRING p-value retriever has started cycle {k} of bottle {nr1}-{nr2}.")

		for j in range(1, (num_cycles+1)):

			list_of_random_stringids = []
			this_line = [j, j]
			for i in range(1, 11):
				this_random_selection = random.choice(list_of_bottle_proteins)
				this_string_id = convert_uniprot[this_random_selection]

				if this_random_selection in name_uniprot:
					this_protein_name = name_uniprot[this_random_selection]
				else:
					this_protein_name = "n/a"

				list_of_random_stringids.append(taxid + "." + this_string_id)
				this_line.append(this_protein_name)
				this_line.append(this_string_id)

			string_api_url = "https://string-db.org/api"
			output_format = "tsv-no-header"
			method = "ppi_enrichment"

			request_url = "/".join([string_api_url, output_format, method])

			params = {
				"identifiers": "%0d".join(list_of_random_stringids),  # list of my selected proteins in STRING correct IDs
				"species": int(taxid),  # species NCBI identifier, ex. human 9606
				"caller_identity": "zdultester"  # my random app name
			}

			## Calling STRING
			response = requests.post(request_url, data=params)

			## Parse and print the respons Parse and print the responsee
			# print("The following IDs have been used: ", ','.join(map(str, list_of_random_stringids)) )

			for line in response.text.strip().split("\n"):
				pvalue = line.split("\t")[5]
				this_line[1] = pvalue
				p_values.append(float(pvalue))

			log_calls.append([nr1, nr2, k, j, pvalue, request_url, response.text.strip()])

			# Print P-value to the console (inactive)
			# print("P-value:", pvalue)

			responses_array.append(this_line)

		pvs = pd.Series(p_values, index=range(1, 11))
		responses_array.append(["SUM", "", "Min (p-val)", pvs.min(), "Max (p-val)", pvs.max(), "Mean (p-val)", pvs.mean()])
		print(f"Parser summary for cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {pvs.min()}.")

		export_filename = "export/pvalues_" + str_goid + "_c" + str(num_cycles).zfill(2) + "_r" + str(num_request_per_cycle).zfill(2) +\
						  "_bottle-" + str(nr1) + "-" + str(nr2) + "_" + taxid + "_" + current_time_abbrev + ".tsv"
		counter = 0

		with open(export_filename, mode='a') as export_file:
			writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

			for line in responses_array:
				counter += 1
				writer.writerow(line)

		# print(f"Parser have {counter} lines wrote in {export_filename} (merged) file.")

# Writing log file
log_filename = "logs/pvalues_" + str_goid + "_c" + str(num_cycles).zfill(2) + "_r" + str(num_request_per_cycle).zfill(2) + "_bottle-" +\
			   str(nr1) + "-" + str(nr2) + "_" + taxid + "_" + current_time_abbrev + ".tsv"

with open(log_filename, mode='a') as log_file:
	writer = csv.writer(log_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

	log_counter = 0
	for line in log_calls:
		writer.writerow(line)
		log_counter += 1

print(f"Parser has written {log_counter}: lines in {log_filename}.")

