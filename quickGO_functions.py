
# QuickGO Parser
#
# What this file do?
# Containing the variables that more than one script use
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

import csv
import random
import time
import logging
from datetime import datetime

# Creating timestamp for output filename
now = datetime.now()
current_time_abbrev = now.strftime("%Y%m%d-%H%M%S-%f")

# LOGGING SETTINGS
dir_log = "logs/"
log_filename = dir_log + "quickgo_" + current_time_abbrev + ".tsv"

# START LOGGING
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename, level=logging.DEBUG)


def TimeNow(str_now, name_of_script=False, start_time=False):
	now = datetime.now()
	time_abb = now.strftime("%Y-%m-%d - %H:%M:%S (%f)")

	if name_of_script:
		if str_now == "start":
			print(f"Runtime of '{name_of_script}' is at {str_now} phase at {time_abb}")
			return True

	elif start_time == False:
		print(f"Runtime is phase {str_now} at {time_abb}")
		return True

	if str_now == "end" or name_of_script == False:

		runtime = time.time() - start_time
		hours = runtime // 3600
		temp = runtime - 3600 * hours
		minutes = temp // 60
		seconds = temp - 60 * minutes

		if name_of_script:
			print(f"Runtime of '{name_of_script}' was", '%d hrs %d mins %d secs' % (hours, minutes, seconds),
				"(" + str(float("{:.5f}".format(runtime))) + ")")
		else:
			print(f"Script running is phase {str_now}, since start", '%d hrs %d mins %d secs elapsed' % (hours, minutes, seconds),
				"(" + str(float("{:.5f}".format(runtime))) + ")")


	return True


def Sleep(func):
	def wrapper():
		sleep_random_time = random.randrange(5, 15)
		print(f"Script is going to sleep {sleep_random_time} seconds.")
		time.sleep(sleep_random_time)
		return func()
	return wrapper


@Sleep
def SleepWakeUp():
	print("Script woke up and continue the process.")


def GOSlimRequestURL(this_go_id, this_taxid):

    go_id_include = this_go_id.replace("_", "%3A")
	# 4932%2C4896%2C9606%2C7227%2C7955%2C6239%2C3702%2C284812%2C559292
    requestURL = "https://" \
                 + "www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" \
                 + "selectedFields=geneProductId&" \
                 + "selectedFields=symbol&" \
                 + "selectedFields=goId&" \
                 + "selectedFields=goName&" \
                 + "selectedFields=taxonId&" \
                 + "selectedFields=synonyms&" \
                 + "selectedFields=goAspect&" \
                 + "goId=" + go_id_include + "&" \
                 + "goUsage=descendants&" \
                 + "goUsageRelationships=is_a%2Cpart_of%2Coccurs_in&" \
                 + "taxonId=" + str(this_taxid) + "&" \
                 + "taxonUsage=descendants"

    return requestURL


def WriteTSVFile(go_id, tax_id, export_filename, write_this_array, split=False, this_mode='a'):

	counter = 0

	with open(export_filename, mode=this_mode) as export_file:
		writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

		for line in write_this_array:
			if split:
				line = line.split(split)

			counter += 1

			try:
				writer.writerow(line)
			except:
				log_this = "GO id: " + go_id + "TaxID: " + str(tax_id) + "Counter: " + counter
				logging.error(log_this)

	log_this = "GO id: " + go_id + "TaxID: " + str(tax_id) + "Filename: " + export_filename + "successful."
	logging.info(log_this)

	return counter
