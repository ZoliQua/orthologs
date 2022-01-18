
# Ortholog Parser / STRING DB Reader
#
# What this file do?
# Containing the variables that more than one script use
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

from datetime import datetime

# Creating timestamp for output filename
now = datetime.now()
current_time_abbrev = now.strftime("%Y%m%d-%H%M%S-%f")

###############################
# Taxon files: list and dicts #
###############################

# Taxon list (id list) tuple
taxon_list = ('9606', '7955', '6239', '3702', '7227', '4896', '4932')

# Taxon dict taxid :: taxon_name
taxon_dict = {	'9606': 'H. sapiens',
				'7955': 'D. rerio',
				'6239': 'C. elegans',
				'3702': 'A. thaliana',
				'7227': 'D. melanogaster',
				'4896': 'S. pombe',
				'4932': 'S. cerevisiae' }

# Taxon column names
taxon_dict_go = {	'9606': 'H. sapiens Hit',
					'7955': 'D. rerio Hit',
					'6239': 'C. elegans Hit',
					'3702': 'A. thaliana Hit',
					'7227': 'D. melanogaster Hit',
					'4896': 'S. pombe Hit',
					'4932': 'S. cerevisiae Hit'
					 }
