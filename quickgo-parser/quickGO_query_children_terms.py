
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

LogAndPrint("Script started.")

for go_name, go_id in GOslim_dict.items():

	# Change GO id "_" => ":"
	go_id_call = go_id.replace("_", ":")
	# Get Children of this GO-id
	go_children = Children(go_id_call)
	if not go_children.successful_run and not go_children.exist and not go_children.no_children:
		LogAndPrint(f"### Error. Repeat once {go_id}. Error ###")
		go_children = Children(go_id_call)
		if not go_children.successful_run:
			LogAndPrint(f"### Error. Second run of {go_id}. Error ###")
		else:
			LogAndPrint(f"### Error 2 OK. Second run of {go_id}. ###")
	if go_children.exist:
		LogAndPrint(f"Skipping. {go_id} file exist.")
	if go_children.no_children:
		LogAndPrint(f"Skipping. {go_id} has no children terms.")

LogAndPrint("Script finished.")
