
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
import signal
import time
import requests  # python -m pip install requests

# maximalise file-read size
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


def ParseGODataFrame(go_dataframe, column_name_taxon, column_name_hm_value, num_proteins):
	
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

	if len(list_of_bottle_proteins) < num_proteins:
		# Print & Log warning
		print(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		logging.warning(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		return False

	else:
		return_array = {"protein_hm_array": protein_hm_array, "list_of_bottle_proteins": list_of_bottle_proteins}
		return return_array

# Custom timeout exception
class TimeoutError(Exception): pass

#  Call this function exceeds timeout
def handler(signum, frame):
    raise TimeoutError()

#  Function timeout decorator
def time_out(interval, doc):
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(interval)       #  Interval seconds to send SIGALRM signals to the process
                result = func(*args, **kwargs)
                signal.alarm(0)              #  After the function is executed after the specified time is executed, close the Alarm alarm clock
                return result
            except TimeoutError as e:
                #  Capture the timeout exception, what to do
                print("The function failed to run due to timeout, func:<%s>" % doc)
                logging.warning("The function failed to run due to timeout, func:<%s>" % doc)
        return wrapper
    return decorator

@time_out(1, "Function call")
def request_in_time(req, par):
    #print("task1 start")
    try:
        response = requests.post(req, data = par)
    except:
        response = False
    return response

