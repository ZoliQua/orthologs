
# Ortholog Parser / STRING DB Reader
#
# What this file do?
# Containing the variables that more than one script use
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

##############################################################
##
## The script prints out the p-value of STRING protein-protein
## interaction enrichment method for the given set of proteins 
##
## Requires requests module:
## type "python -m pip install requests" in command line (win)
## or terminal (mac/linux) to install the module
## from STRING API documentation
##
## URL: https://string-db.org/help/api/
##
##############################################################

import requests

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "ppi_enrichment"

##########################
## Construct the request #
##########################

request_url = "/".join([string_api_url, output_format, method])

###################
## Set parameters #
###################

# Test set of genes in human (taxid: 9606)
# my_taxid = 9606
# my_genes = ['9606.ENSP00000004531', '9606.ENSP00000231706', '9606.ENSP00000291900', '9606.ENSP00000294353',
#             '9606.ENSP00000346693', '9606.ENSP00000358831', '9606.ENSP00000359956', '9606.ENSP00000360583',
#             '9606.ENSP00000363417', '9606.ENSP00000370128', '9606.ENSP00000372390', '9606.ENSP00000424123',
#             '9606.ENSP00000477602']

# Test set of genes in fly (taxid: 7227)
my_taxid = 7227
my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
            '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']

params = {
    "identifiers":      "%0d".join(my_genes),   # your proteins
    "species":          my_taxid,                   # species NCBI identifier
    "caller_identity":  "tester_zdul"           # your app name
}

###################
## Call STRING DB #
###################

response = requests.post(request_url, data=params)

##
## Parse and print the response (incl. p-value) into the console
##

#################################
## Output fields of STRING API ##
#################################
# Field 							Description
#
# number_of_nodes				number of proteins in your network
# number_of_edges 				number of edges in your network
# average_node_degree			mean degree of the node in your network
# local_clustering_coefficient	average local clustering coefficient
# expected_number_of_edges		expected number of edges based on the nodes degrees
# p_value						significance of your network having more interactions than expected

for line in response.text.strip().split("\n"):
    number_of_nodes = line.split("\t")[0]
    number_of_edges = line.split("\t")[1]
    average_node_degree = line.split("\t")[2]
    local_clustering_coefficient = line.split("\t")[3]
    expected_number_of_edges = line.split("\t")[4]
    pvalue = line.split("\t")[5]

    print("")
    print(f"The query of the following genes ({len(my_genes)}):")
    print("="*40)
    print(", ".join(my_genes))
    print("="*40)
    print("Number_of_nodes".ljust(30), number_of_nodes)
    print("Number_of_edges".ljust(30), number_of_edges)
    print("Average_node_degree".ljust(30), average_node_degree)
    print("Local_clustering_coefficient".ljust(30), local_clustering_coefficient)
    print("Expected_number_of_edges".ljust(30), expected_number_of_edges)
    print("P-value:".ljust(30), pvalue)
