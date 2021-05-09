
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
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename, level=logging.INFO)


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


call_counter = 0


def Sleep(func):
	"""Sleeps a given time and returns"""
	global call_counter

	def wrapper():
		sleep_random_time = random.randrange(5, 15)
		LogAndPrint(f"Script is going to sleep {sleep_random_time} seconds.", False)
		time.sleep(sleep_random_time)
		return func(1)

	def counter():
		global call_counter
		call_counter += 1
		return func(0)

	if call_counter >= 100:
		call_counter = 0
		return wrapper
	else:
		return counter


@Sleep
def SleepWakeUp(state=0):
	"""Returns from sleep"""
	if state == 1:
		LogAndPrint("Script woke up and continue the process.", False)


def LogAndPrint(text, is_printing=True, level="info"):
	"""Takes a text, that logs into a log file and print into the console"""
	if level == "info":
		logging.info(text)
	else:
		logging.debug(text)
	if is_printing:
		print(text)
	return True


def GOSlimRequestURL(this_go_id, this_taxid):

	# Change "_" into ":" in the GO term name
    go_id_include = this_go_id.replace("_", "%3A")

	# All seven taxids to include into taxonID place
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


def GetRequest(goid):
	"""Takes a request for QuickGO, returns page content"""
	url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + goid + "/children"
	hdr = {'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8', 'User-Agent': "Magic Browser"}
	page = rqs.Request(url, headers=hdr)

	return rqs.urlopen(page).read()


list_of_goterms = []
dict_of_goterms = {}
export_to_tsv = []


def GetChildren(goid, level, parent_goid=None):

	SleepWakeUp()

	content = GetRequest(goid)
	json_content = json.loads(content)

	if level == 0:
		LogAndPrint(f'Start {goid} request from EBI.')
		parent_goid = goid

	for children in json_content['results'][0]['children']:
		if children['hasChildren']:
			if children['id'] not in list_of_goterms:
				LogAndPrint("\t"*level
				            + "(lv-" + str(level) + ") "
				            + 'Children of ' + children['id'] + " (" + children['name'] + "):")
				GetChildren(children['id'], level+1, children['id'])
				list_of_goterms.append(children['id'])
				dict_of_goterms[children['id']] = children['name']
				line = [children['id'], children['name'], children['relation'], 1, parent_goid]
				export_to_tsv.append(line)
			else:
				LogAndPrint("\t"*level + "(lv-" + str(level) + ") " + children['id'] + ' (children above)')
		else:
			if children['id'] not in list_of_goterms:
				list_of_goterms.append(children['id'])
				dict_of_goterms[children['id']] = children['name']
				line = [children['id'], children['name'], children['relation'], 0, parent_goid]
				export_to_tsv.append(line)
				LogAndPrint("\t"*level + "(lv-" + str(level) + ") " + children['id'] + ": " + children['name'])
			else:
				LogAndPrint("\t"*level + "(lv-" + str(level) + ") " + children['id'] + " (repeat)")

	return True


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


def WriteTSVFileChildren(export_filename, write_this_array, this_mode='a'):

	counter = 0

	with open(export_filename, mode=this_mode) as export_file:
		writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

		for line in write_this_array:
			counter += 1
			writer.writerow(line)

	log_this = "Filename: " + export_filename + " successful."
	LogAndPrint(log_this)

	return counter
