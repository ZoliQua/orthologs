
import sys
import os
from PIL import Image

image_folder = sys.argv[1]
output_folder = sys.argv[2]

print(f'Target directory: {image_folder}')
print(f'Destination directory: {output_folder}')

if not os.path.exists(output_folder):
	os.makedirs(output_folder)
	print(f'A folder named {output_folder} has been created.')
else:
	print(f'A folder named {output_folder} already exists.')

images = os.listdir(image_folder)

# print(images)

for image in images:
	img = Image.open(f'{image_folder}/{image}')
	cleaned_name = os.path.splitext(image)[0]
	img.save(f'{output_folder}/{cleaned_name}.png', 'png')
	print(f'Script converted {image} to {output_folder}/{cleaned_name}.png')

print("Run was successful.")