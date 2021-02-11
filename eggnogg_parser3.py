#This file reads the eggNOG db to parse it to the 7 speceies

import csv
import sys

csv.field_size_limit(sys.maxsize)

filename = "data/eggnog/2759_members.tsv"
taxon_dict = {'9606': 'H. sapiens', '7955': 'D. rerio', '6239': 'C. elegans', '3702': 'A. thaliana', '7227': 'D. melanogaster', '4896': 'S. pombe', '284812': 'S. pombe', '559292': 'S. cerevisiae'}
uniprot_list = {'9606': [], '7955': [], '6239': [], '3702': [], '7227': [], '4896': [], '4932': []}
eggnog_taxlist = ["9606", "7955", "6239", "3702", "7227", "4896", "4932"]
final_list = {}
export_text = ""

# for i in eggnog_taxlist:
#     print(i)

with open(filename, newline='') as f:
    reader = csv.DictReader(f, fieldnames= ('eggnogmain', 'groupid','numprot','numspec','prots','specs' ), delimiter='\t')
    counter = 0
    try:
        for row in reader:

            row_counter = 0

            this_prots = row['prots'].split(",")
            this_specs = row['specs'].split(",")

            for i in this_specs:
                if i in eggnog_taxlist:
                    row_counter += 1

            if row_counter == 0:
                continue
            else:
                final_list[row['groupid']] = {'9606': [], '7955': [], '6239': [], '3702': [], '7227': [], '4896': [], '4932': []}
                for i in this_prots:
                    splited = i.split(".")

                    if splited[0] in eggnog_taxlist:
                        final_list[row['groupid']][splited[0]].append(i[5:])
                        uniprot_list[splited[0]].append(i[5:])


    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

counter = 0

with open('data/eggnogdb_7_species_by_groups.tsv', mode='w') as export_file:
    writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    for groupid in final_list:

        counter += 1
        this_line = [groupid]

        countspec = 0

        this_line.append(countspec)

        for i in eggnog_taxlist:
            this_line.append(",".join(final_list[groupid][i]))
            if len(final_list[groupid][i]) > 0:
                countspec += 1

        this_line[1] = countspec

        writer.writerow(this_line)

        # if counter == 50:
        #     break

for i in eggnog_taxlist:

    filename = 'data/eggnogdb_'+i+'_protein_list.tsv'
    print(filename)

    with open(filename, mode='w') as export_file:
        writer = csv.writer(export_file, delimiter='\r', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # this_line = []
        # this_line.append(uniprot_list[i])
        writer.writerow(uniprot_list[i])

# print(uniprot_list)


#print(uniprot_list)
