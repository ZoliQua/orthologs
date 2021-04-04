##############################################################
## The script prints out the p-value of STRING protein-protein
## interaction enrichment method for the given set of proteins 
##
## Requires requests module:
## type "python -m pip install requests" in command line (win)
## or terminal (mac/linux) to install the module
## from STRING API documentation
##############################################################

import requests

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "ppi_enrichment"

##
## Construct the request
##

request_url = "/".join([string_api_url, output_format, method])

##
## Set parameters
##

# my_genes = ['9606.ENSP00000216911', '9606.ENSP00000237654', '9606.ENSP00000256442', '9606.ENSP00000261819', '9606.ENSP00000282572', '9606.ENSP00000287598',
#             # '9606.ENSP00000288207', '9606.ENSP00000302530', '9606.ENSP00000302898', '9606.ENSP00000313950', '9606.ENSP00000314004', '9606.ENSP00000315743',
#             # '9606.ENSP00000339109', '9606.ENSP00000344635', '9606.ENSP00000357858', '9606.ENSP00000365210', '9606.ENSP00000394394', '9606.ENSP00000396755',
#             '9606.ENSP00000426654', '9606.ENSP0000047825']

# Test set of genes in human (taxid: 9606)
my_genes = ['9606.ENSP00000004531', '9606.ENSP00000231706', '9606.ENSP00000291900', '9606.ENSP00000294353', '9606.ENSP00000346693', '9606.ENSP00000358831',
            '9606.ENSP00000359956', '9606.ENSP00000360583', '9606.ENSP00000363417', '9606.ENSP00000370128', '9606.ENSP00000372390', '9606.ENSP00000424123',
            '9606.ENSP00000477602']

# Test set of genes in fly (taxid: 7227)
# my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
#             '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']

params = {
    "identifiers" : "%0d".join(my_genes), # your proteins
    "species" : 9606, # species NCBI identifier
    "caller_identity" : "testerrr" # your app name
}

##
## Call STRING DB
##

response = requests.post(request_url, data=params)

##
## Parse and print the response (incl. p-value) into the console
##

for line in response.text.strip().split("\n"):
    pvalue = line.split("\t")[5]
    print("P-value:", pvalue)