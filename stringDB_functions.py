
# File for funcitions
# Import libraries
import csv
import sys

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
