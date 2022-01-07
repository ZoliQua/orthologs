import pandas as pd
import os
import re
from bs4 import BeautifulSoup

# Read the TSV file
data_file = 'go-0007049-ordered-2022-01.tsv'
data = pd.read_csv(data_file, sep='\t')

# Count the overlapping groups
groups = {}
for index, row in data.iterrows():
    key = ''
    for col in range(14, 21):  # Columns from 15 to 21 (index 14 to 20)
        if not pd.isna(row[col]):
            key += chr(65 + col - 14)  # Convert column number to corresponding uppercase letter (A to G)
    if key:
        if key not in groups:
            groups[key] = 0
        groups[key] += 1

# Select the appropriate SVG file based on the number of shapes
shapes = len(data.columns[14:21])
svg_file = f'source/ortholog_venn_{shapes}.svg'

# Read the SVG file
with open(svg_file, 'r') as f:
    svg_data = f.read()

# Modify the SVG file using BeautifulSoup
soup = BeautifulSoup(svg_data, 'html.parser')
texts = soup.find_all('text')

# Update the overlapping counts in the SVG file
for text in texts:
    if text.string in groups:
        count = groups[text.string]
        text.string = str(count)
    elif re.match(r'^[A-G]+$', text.string):  # Check if the text contains only uppercase letters from A to G
        text.string = '0'

    # Update the header names
    if text.string.startswith('Name'):
        col = ord(text.string[-1]) - 65  # Convert the last character to a column index
        header_name = data.columns[14 + col]
        text.string = header_name

# Write the modified SVG file
output_file = f'output/ortholog_venn_modified_{shapes}_hit.svg'
with open(output_file, 'w') as f:
    f.write(str(soup))