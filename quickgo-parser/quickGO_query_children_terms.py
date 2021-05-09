
# QuickGO Parser
#
# What this file do?
# Containing the variables that more than one script use
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

from quickGO_parse_GOslim import *
from quickGO_functions import *

for go_name, go_id in GOslim_dict.items():

	GetChildren(go_id, 0)
	filename = "export/" + this_go.replace(":", "_") + ".tsv"
	WriteTSVFileChildren(filename, export_to_tsv, "w")

LogAndPrint("Finished")
