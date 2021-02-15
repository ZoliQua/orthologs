
# UNIPROT PARSER v1.0
#
# What this file do?
# This file converts downloaded id_mapping files (from data/uniprot folder) to a reduced file size. Filtering out STRING and eggNOG db related lines.
#
# Source folder for files (data/unirprot): https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
# Source files downloaded on 11/02/2021
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import csv
import sys

csv.field_size_limit(sys.maxsize)

filenames = ["data/uniprot/ARATH_3702_idmapping.dat", "data/uniprot/CAEEL_6239_idmapping.dat", "data/uniprot/DANRE_7955_idmapping.dat", "data/uniprot/DROME_7227_idmapping.dat", "data/uniprot/HUMAN_9606_idmapping.dat", "data/uniprot/SCHPO_284812_idmapping.dat", "data/uniprot/YEAST_559292_idmapping.dat"]
taxon_dict = {'9606': 'H. sapiens', '7955': 'D. rerio', '6239': 'C. elegans', '3702': 'A. thaliana', '7227': 'D. melanogaster', '4896': 'S. pombe', '284812': 'S. pombe', '559292': 'S. cerevisiae'}
uniprot_list = {'9606': [], '7955': [], '6239': [], '3702': [], '7227': [], '4896': [], '4932': []}
eggnog_taxlist = ["9606", "7955", "6239", "3702", "7227", "4896", "4932"]
write_lines_all = []

for this_file in filenames:

    counter = 0
    write_lines = []
    filename_taxid_split = this_file.split("_")
    this_taxid = filename_taxid_split[1]

    with open(this_file, newline='') as f:
        reader = csv.DictReader(f, fieldnames= ( 'uniprot', 'db','second' ), delimiter='\t')
        counter = 0
        try:
            for row in reader:
                counter += 1
                # Filtering eggNOG DB out
                if row['db'] == 'eggNOG':
                    write_lines.append("\t".join([row['uniprot'], row['db'], this_taxid, row['second']]))
                    write_lines_all.append("\t".join([row['uniprot'], row['db'], this_taxid, row['second']]))
                # Filtering STRING DB out
                if row['db'] == 'STRING':
                    write_lines.append("\t".join([row['uniprot'], "convert", this_taxid, row['second'][5:]]))
                    write_lines_all.append("\t".join([row['uniprot'], "convert", this_taxid, row['second'][5:]]))

                # More filters can be added easily by repeating the condion above with different terms.

                # if counter == 10000:
                #     break

        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

        print("Parser have", counter, "lines processed.")

    counter = 0

    export_filename = "data/uniprot_convert_"+this_taxid+".tsv"

    with open(export_filename, mode='w') as export_file:
        writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting = csv.QUOTE_MINIMAL)

        for line in write_lines:

            counter += 1
            this_line = line.split("\t")
            writer.writerow(this_line)

        print("Parser have ", counter, " lines wrote in", export_filename, "file.")

export_filename_merged = "data/uniprot_convert_merged.tsv"

with open(export_filename_merged, mode='w') as export_file:
    writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting = csv.QUOTE_MINIMAL)

    for line in write_lines_all:

        counter += 1
        this_line = line.split("\t")
        writer.writerow(this_line)

    print("Parser have ", counter, " lines wrote in", export_filename, "(merged) file.")