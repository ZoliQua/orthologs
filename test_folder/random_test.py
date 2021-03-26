#
# In this file I tested random generators and plotting in matplotlib.
#

import random
import string
from statistics import mean
from datetime import datetime

#
# Importing MatPlotLib
#

from matplotlib import pyplot as plt
from matplotlib import image as mpimage

#
# Here comes the code
#

# List of sources for random selection
list_of_sources = [string.ascii_letters,
				   string.ascii_lowercase,
				   string.ascii_uppercase,
				   string.digits,
				   string.hexdigits,
				   string.octdigits,
				   string.punctuation,
				   string.printable]

# Randomly selecting the source
selected_source = random.choice(list_of_sources)
print(f"Selecting from the following letters \"{selected_source}\".")

# Randomly selecting the a letter from the source
selected_random_letter = random.choice(selected_source)
print(f"Random function selected the following letter: \"{selected_random_letter}\" at position {selected_source.find(selected_random_letter)}.")

#
# Creating a random 10 digit long int
#

random_int = ""

# Creating the 10 digit long int with for loop
for i in range(1, 11):
	random_int += str(random.randrange(1, 10))

# Print the 10 digit long int
# print(random_int)

#
# Testing random function
# Iterating 9 times a for loop that do the following:
# 1.) Run a for loop 100 times to run a random number selector
# 2.) Random number selector each time run 10 000 (10k) times:
# 		randomly selecting a number between 1 and 9
# 		then count it in an array (int_list)
# 		after 10k: select
# 					the min count number,
# 					max count number,
# 					calculate the max/min ratio
# 3.) Collect the max/min ratio in a list
# 4.) Plot with mathplotlib the array of the max/min ratio
# 5.) Save the output file in dir: jpgImageDir
# In summary this script creates 9 different JPG files in the output folder.
#


for j in range(1, 10):

	list_of_min = []
	list_of_max = []
	list_of_minmax = []

	for i in range(1, 100):

		# Creating initial list with 10 elements
		int_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		# Iterating 10 000 times the random selection from 1-9, then counting them into an array.
		for i in range(1, 10_000):
			int_list[random.randrange(1, 10)] += 1

		# We don't need number 0 element. It has a value of 0 we know.
		int_list.pop(0)

		# Print out the values to see
		# print(int_list)
		# print(f"Minimum number is: {min(int_list)}.")
		# print(f"Maximum number is: {max(int_list)}.")

		list_of_precents = []
		minmax = float("{:.6f}".format(max(int_list) / min(int_list)))

		# Iterating the current list to determine the min, the max and max/min ratio values
		for int_num, count in enumerate(int_list):
			int_num += 1
			percent = float("{:.6f}".format(count / min(int_list)))
			list_of_precents.append(percent)
			# print(f"{int_num}: {count}, {percent} ")

		# Writing to arrays the min, the max and max/min ratios
		list_of_min.append(min(list_of_precents))
		list_of_max.append(max(list_of_precents))
		list_of_minmax.append(minmax)

	#
	# Plotting with MatPlotLib
	#

	# Creating timestamp for output plot filename
	now = datetime.now()
	current_time_abbrev = now.strftime("%Y%m%d-%H%M%S-%f")

	# Text for the x axis
	plt.xlabel("Number of iterations")
	# Text for the y axis
	plt.ylabel("Max/Min ratio")
	# Title of the plot
	plt.title("Random numbers")
	# Plotting
	plt.plot(list_of_minmax, "-", label="test", linewidth=3, color="orange")
	# Save the file
	plt.savefig('jpgImageDir/line_plot_' + current_time_abbrev + '.jpg', dpi=300)
	# Show the file
	plt.show()
	# Print the date into the console
	print(list_of_minmax)







