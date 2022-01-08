import pandas as pd
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

# Read the TSV file
file_path = "go-0007049-ordered-2022-01.tsv"
df = pd.read_csv(file_path, sep='\t')

# Initialize the count dictionary
overlap_counts = defaultdict(int)

# Iterate through the rows and count the overlaps
for index, row in df.iterrows():
    count = 0
    for i in range(8, 15):
        if pd.notna(row[i]):
            count += 1
    overlap_counts[count] += 1

# Function to modify the SVG file
def modify_svg_file(input_file, output_file, overlap_counts):
    tree = ET.parse(input_file)
    root = tree.getroot()

    # Update the SVG elements with overlap counts
    for key, value in overlap_counts.items():
        if key >= 2:
            element_id = ''.join(chr(ord('A') + i) for i in range(key))
            for element in root.findall(".//{http://www.w3.org/2000/svg}text"):
                if element.text == element_id:
                    element.text = str(value)
                    break

    # Save the modified SVG file
    tree.write(output_file)

# Determine the appropriate SVG file to modify
num_shapes = max(overlap_counts.keys())
source_svg_file = f"source/ortholog_venn_{num_shapes}.svg"
output_svg_file = f"output/ortholog_venn_{num_shapes}_modified.svg"
Path("output").mkdir(exist_ok=True)

# Modify the SVG file with the overlap counts
modify_svg_file(source_svg_file, output_svg_file, overlap_counts)
