
# EGGNOG RANDOM SELECTOR
#
# What this file do?
# x
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import os
import csv
import sys
import random

csv.field_size_limit(sys.maxsize)

taxid = "9606"

filename = "data/uniprot_convert_merged.tsv"
write_lines_all = {}
convert_uniprot = {}

with open(filename, newline='') as f:
    reader = csv.DictReader(f, fieldnames= ('uniprot', 'db', 'taxid', 'selector'), delimiter='\t')
    counter = 0
    try:
        for row in reader:
            counter += 1

            # Filtering STRING convert out
            if row['db'] == 'convert' and :

                if row['uniprot'] in convert_uniprot:
                    continue
                else:
                    convert_uniprot[row['uniprot']] = row['selector']

            # if counter == 10000:
            #      break

    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

    print("Parser have", counter, "lines processed from eggNOG merged file.")

    random_entry = random.choice(list(convert_uniprot.values()))
    print(random_entry)

