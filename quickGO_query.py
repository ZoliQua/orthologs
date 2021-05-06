
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

# Print start time to the console
start_time = time.time()
this_file = "quickGO_query.py"
TimeNow("start", this_file)


folder = "data/go/"

for go_name, go_id in GOslim_dict.items():

    phase_name = "\'GO Request for " + go_id.replace("_", ":") + " (" + go_name + ")\'"
    TimeNow(phase_name)

    this_filename = folder + go_id + ".tsv"

    r = requests.get(GOSlimRequestURL(go_id), headers={"Accept": "text/tsv"})

    if not r.ok:
        r.raise_for_status()
        # sys.exit()
        continue

    responseBody = r.text
    responseBodyArray = responseBody.split("\n")

    phase_name = "\'GO Process " + go_id.replace("_", ":") + " (" + go_name + "), lines " \
                 + len(responseBodyArray)-2 + "\'"
    TimeNow(phase_name, False, start_time)

    WriteTSVFile(this_filename, responseBodyArray, "\t")

    # SleepWakeUp()
    # break

# Print end time to the console
TimeNow("end", this_file, start_time)
