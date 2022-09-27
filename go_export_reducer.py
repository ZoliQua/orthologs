# GO PARSER v1.0
#
# What this file do?
# This file converts downloaded GO export file (from data/go folder) to a reduced file size.
# Filtering out tax_id, sub_GO term, and uniprot ids
# Then print it to the console.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import csv

filename = "data/go/go_reg_of_cc.tsv"
taxon_list = []
with open(filename, newline='') as f:
    reader = csv.DictReader(f, fieldnames= ('type', 'uniprot', 'name', 'subgo', 'longname', 'taxid'), delimiter='\t')
    counter = 0
    try:
        for row in reader:
            if row['type'] != 'UniProtKB':
                continue
            print(row['taxid'], row['subgo'], row['uniprot'])
            counter += 1
            if row['taxid'] in taxon_list:
                continue
            else:
                taxon_list.append(row['taxid'])
            # if counter == 100:
            #    break
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

print(taxon_list)