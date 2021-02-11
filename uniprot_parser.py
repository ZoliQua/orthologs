
#This file converts uniprod downloaded mapping files to a reduced file size

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

    with open(this_file, newline='') as f:
        reader = csv.DictReader(f, fieldnames= ( 'uniprot', 'db','second' ), delimiter='\t')
        counter = 0
        try:
            for row in reader:
                counter += 1
                if row['db'] == 'eggNOG':
                    write_lines.append("\t".join([row['uniprot'], row['db'], row['second']]))
                    write_lines_all.append("\t".join([row['uniprot'], row['db'], row['second']]))

                if row['db'] == 'STRING':
                    write_lines.append("\t".join([row['uniprot'], "convert", row['second'][5:]]))
                    write_lines_all.append("\t".join([row['uniprot'], "convert", row['second'][5:]]))

                # if counter == 10000:
                #     break

        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

        print("Parser have", counter, "lines processed.")

    counter = 0

    export_filename_taxid = this_file.split("_")
    export_filename = "data/uniprot_convert_"+export_filename_taxid[1]+".tsv"

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