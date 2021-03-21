
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
import requests ## python -m pip install requests

csv.field_size_limit(sys.maxsize)

taxid = "9606"
human_ids = ["Q5TKA1", "Q9H1A4", "Q9UJX5", "Q9UJX4", "Q9UJX3", "Q13535", "O14965", "O60566", "O43684", "P14635", "P24863", "P24385", "P51946", "Q9UK58", "O75909", "P30260", "Q13042", "Q9UJX2", "Q99741", "Q8IZL9", "Q9BWU1", "Q9H211", "Q96JB5", "P61024", "P68400", "Q14186", "Q9Y6D9", "Q13257", "Q6NUQ1"]

filename = "data/uniprot_convert_merged.tsv"
write_lines_all = {}
convert_uniprot = {}

with open(filename, newline='') as f:
    reader = csv.DictReader(f, fieldnames= ('uniprot', 'db', 'taxid', 'selector'), delimiter='\t')
    counter = 0
    try:
        for row in reader:

            # Filtering STRING convert out
            if row['db'] == 'convert' and row['taxid'] == taxid:

                if row['uniprot'] in convert_uniprot:
                    continue
                else:
                    convert_uniprot[row['uniprot']] = row['selector']
                    counter += 1

            # if counter == 10000:
            #      break

    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

    print("Parser have", counter, "lines processed from eggNOG merged file.")

for j in range(1,11):

    list_of_random_stringids = []
    for i in range(1,8):
        list_of_random_stringids.append('9606.'+convert_uniprot[random.choice(list(human_ids))])

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "ppi_enrichment"

    request_url = "/".join([string_api_url, output_format, method])

    params = {
        "identifiers" : "%0d".join(list_of_random_stringids),   # list of my selected proteins in STRING correct IDs
        "species" : int(taxid),                                 # species NCBI identifier, ex. human 9606
        "caller_identity" : "zdultester"                          # my random app name
    }

    ## Calling STRING

    response = requests.post(request_url, data=params)

    ## Parse and print the respons Parse and print the responsee

    print("The following IDs have been used: ", ','.join(map(str, list_of_random_stringids)) )

    for line in response.text.strip().split("\n"):
        pvalue = line.split("\t")[5]
        print("P-value:", pvalue)
