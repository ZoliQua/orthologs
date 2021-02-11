
import csv

filename = "data/go_reg_of_cc.tsv"
taxon_dict = {'9606': 'H. sapiens', '7955': 'D. rerio', '6239': 'C. elegans', '3702': 'A. thaliana', '7227': 'D. melanogaster', '4896': 'S. pombe', '284812': 'S. pombe', '559292': 'S. cerevisiae'}
uniprot_list = {'9606': [], '7955': [], '6239': [], '3702': [], '7227': [], '4896': [], '284812': [], '559292': []}
list_uniprot2taxid = {}
list_uniprotall = []
list_uniprotselected = []
list_unip2eggnog = {}
list_eggnog2unip = {}

with open(filename, newline='') as f:
    reader = csv.DictReader(f, fieldnames= ('type', 'uniprot','name','subgo','longname','taxid' ), delimiter='\t')
    counter = 0
    try:
        for row in reader:
            if row['type'] != 'UniProtKB':
                continue
            # print(taxon_dict[row['taxid']], row['subgo'], row['uniprot'])

            counter += 1

            list_uniprot2taxid[row['uniprot']] = row['taxid']
            list_uniprotall.append(row['uniprot'])

            if (row['uniprot'] in uniprot_list[row['taxid']]):
                continue
            else:
                uniprot_list[row['taxid']].append(row['uniprot'])
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))


filename = "data/latest.Eukaryota.tsv"
with open(filename, newline='') as f:
    reader = csv.DictReader(f, fieldnames= ('uniprot', 'eggnog'), delimiter='\t')
    counter = 0
    try:
        for row in reader:

            counter += 1

            # if row['uniprot'] in list_uniprotall:
            #     list_uniprotselected.append(row['uniprot'])
            # else:
            #     continue

            list_unip2eggnog[row['uniprot']] = row['eggnog']
            try:
                list_eggnog2unip[row['eggnog']].append(row['uniprot'])
            except:
                list_eggnog2unip[row['eggnog']] = [row['uniprot']]

            if counter % 200000 == 0:
                print(counter)
                # break

    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

selected_eggnogs = []

for i in list_uniprotall:
    try:
        # print(i, list_uniprot2taxid[i], list_unip2eggnog[i])
        if list_unip2eggnog[i] in selected_eggnogs:
            continue
        else:
            selected_eggnogs.append(list_unip2eggnog[i])
    except:
        continue

#for i in selected_eggnogs:


print(len(selected_eggnogs))

