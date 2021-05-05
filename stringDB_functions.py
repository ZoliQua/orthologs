
# Ortholog Parser / STRING DB Reader
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins in a given bottle.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#
# File for funcitions
# Import libraries
import csv
import sys
import random
import logging

# Maximalize file-read size
csv.field_size_limit(sys.maxsize)

uniprot_2_stringid = {}
uniprot_2_protname = {}
list_of_uniprotids = []

#############
# FUNCTIONS #
#############

#####################
# Read UniProt file #
#####################


###################################################
# Read File Function ##############################
###################################################
#  reads the source tsv file for the given taxid ##
#  takes an str = taxid as an input ###############
#  writes 3 global variables ######################
#  returns a count of read lines ##################
###################################################


def ReadUniprotConvert(taxid):

	global uniprot_2_protname
	global uniprot_2_stringid
	global list_of_uniprotids

	filename = "data/uniprot_convert_" + taxid + ".tsv"

	with open(filename, newline='') as f:
		reader = csv.DictReader(f, fieldnames=('uniprot', 'db', 'taxid', 'selector'), delimiter='\t')
		counter = 0
		try:
			for row in reader:

				# Filtering STRING to convert Uniprot to STRING db id
				if row['db'] == 'convert':
					if row['uniprot'] in uniprot_2_stringid:
						continue
					else:
						uniprot_2_stringid[row['uniprot']] = row['selector']
						list_of_uniprotids.append(row['uniprot'])
						counter += 1

				# Filtering STRING convert out
				if row['db'] == 'Gene_Name':
					if row['uniprot'] in uniprot_2_protname:
						continue
					else:
						uniprot_2_protname[row['uniprot']] = row['selector']

		except csv.Error as e:
			sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

		return counter


#################################
# Write Export Function #########
#################################
#  creates a tsv separated file #
#  takes an array as an input ###
#  returns the count of lines ###
#################################


def WriteExportFile(export_filename, write_this_array):

	counter = 0

	with open(export_filename, mode='a') as export_file:
		writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

		for line in write_this_array:
			counter += 1
			writer.writerow(line)

	return counter


def ParseGODataFrame(go_dataframe, column_name_taxon, column_name_hm_value):
	
	# Protein list from the filtered pandas dataset of this species
	list_of_bottle_proteins = []
	protein_hm_array = {}
	
	for index, row_data in go_dataframe.iterrows():

		protein = str(row_data[column_name_taxon])
		hm_value = row_data[column_name_hm_value]

		# When there is no hit protein for this species: SKIP
		if protein == "nan":
			continue
		# When there are more than one hit protein, randomly selection one
		if protein.find(',') != -1:
			proteins_array = protein.split(",")
			protein_selected = random.choice(proteins_array)
		else:
			protein_selected = protein

		# Adding proteins to a list
		list_of_bottle_proteins.append(protein_selected)
		# Adding poritnes to H/M array
		protein_hm_array[protein_selected] = hm_value

	if len(list_of_bottle_proteins) < 10:
		# Print & Log warning
		print(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		logging.warning(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		return False

	else:
		return_array = {"protein_hm_array": protein_hm_array, "list_of_bottle_proteins": list_of_bottle_proteins}
		return return_array
