#
# QuickGO Annotations shower
#
# What this file does?
# Shows the annoations of a given GO term.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

from quickGO_get_GOslim import *
from quickGO_functions_container import *

LogAndPrint("Script started.")

for go_name, go_id in GOslim_dict.items():

	# Change GO id "_" => ":"
	go_id_call = go_id.replace("_", ":")

	# Get the annotations of this GO ID
	go_children = Annotations(go_id_call)

	break

LogAndPrint("Script finished.")
