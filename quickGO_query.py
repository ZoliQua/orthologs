
# QuickGO Parser
#
# What this file do?
# Containing the variables that more than one script use
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

# taxon ids: 9606, 7955, 6239, 3702, 7227, 4896, 4932, 284812, 559292

from quickGO_parse_GOslim import *
from quickGO_functions import *
import requests, sys
import os

# Print start time to the console
start_time = time.time()
this_file = "quickGO_query.py"
TimeNow("start", this_file)

# List of TaxIDs to check
list_of_taxids = [9606, 7955, 6239, 3702, 7227, 4896, 4932, 284812, 559292]
# Directory of export GO files
dir_export = "data/go/"

for go_name, go_id in GOslim_dict.items():

    # go_id = "GO_0005615"

    # Export filename
    this_filename = dir_export + go_id + ".tsv"

    # Check file status
    if os.path.exists(this_filename):
        os.remove(this_filename)
        logging.info(this_filename + " has been removed.")

    for taxid in list_of_taxids:

        phase_name = "\'GO Request " + str(taxid) + " for " + go_id.replace("_", ":") + " (" + go_name + ")\'"
        logging.info(phase_name + " start.")

        r = requests.get(GOSlimRequestURL(go_id, taxid), headers={"Accept": "text/tsv"})

        if not r.ok:
            logging.error(r)
            r.raise_for_status()
            # sys.exit()
            continue

        responseBody = r.text
        responseBodyArray = responseBody.split("\n")

        # Remove the last line, while it is an empty line
        responseBodyArray.pop()

        # If this is not the first taxid, remove the first line
        if taxid != 9606:
            responseBodyArray.pop(0)

        phase_name = "\'GO Process " + go_id.replace("_", ":") + " - taxid: " + str(taxid) + ", lines " \
                     + str(len(responseBodyArray)) + "\'"

        TimeNow(phase_name, False, start_time)
        logging.info(phase_name + " - end.")

        WriteTSVFile(go_id, taxid, this_filename, responseBodyArray, "\t")

    #   SleepWakeUp()
    # break

# Print end time to the console
TimeNow("end", this_file, start_time)
