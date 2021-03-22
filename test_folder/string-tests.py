
#
# See how many characters are in a given string
#

test_text = "hello!"
print(len(test_text))


#
# Create a list from a given variable string
#

test_listfromtext = list('hello!')
print(test_listfromtext)

#
# Concatenate two different strings
#

test_string1 = "hello"
test_string2 = " "
test_string3 = "world"
test_merged_string = test_string1 + test_string2 + test_string3
print(test_merged_string)

#
# Get copies of a string, then cut the last two characters with substring
# Then try different functions as the names.
#

string = 'hello'
longstring = ('a-' + string + '-b, ') * 5
print(longstring[0:-2])
print("Lower: ", string.lower())
print("Upper: ", string.upper())
print("Capitalize: ", test_merged_string.capitalize())
print("Title: ", test_merged_string.title())

#
# Split a word
#

splitted_text = string.split(sep = 'e')
print("\"", string, "\" splitted as", splitted_text)
splitted_text2 = test_merged_string.split(sep = ' ')
print(splitted_text2)

#
# Substrings and Revers the string
#

print("1:", string[1:])
print("2:", string[2:])
print("3:", string[3:])
print("4:", string[4:])
print(":1", string[:1])
print(":2", string[:2])
print(":3", string[:3])
print(":4", string[:4])
print(":5", string[:5])
print("1:4", string[1:4])
print("2:4", string[2:4])
print("1:5", string[1:5])
print("2:5", string[2:5])
print("2:5:1", string[2:5:1])
print("2:5:2", string[2:5:2])
print("2::2 merged:", test_merged_string[2::2])
print("2:10:2 merged:", test_merged_string[2:10:2])
print("::-1 merged:", test_merged_string[::-1])
print("1::-1 merged:", test_merged_string[1::-1])
print("2::-1 merged:", test_merged_string[2::-1])
print("3::-1 merged:", test_merged_string[3::-1])
print("4::-1 merged:", test_merged_string[4::-1])
print("5::-1 merged:", test_merged_string[5::-1])
print("::-1 merged:", test_merged_string[::-1])
print("2::-1 merged:", test_merged_string[2::-1])
print("1:0:-1 merged:", test_merged_string[1:0:-1])
print("2:0:-1 merged:", test_merged_string[2:0:-1])
print("3:0:-1 merged:", test_merged_string[3:0:-1])
print("4:0:-1 merged:", test_merged_string[4:0:-1])
print("5:0:-1 merged:", test_merged_string[5:0:-1])
print("1:1:-1 merged:", test_merged_string[1:1:-1])
print("2:1:-1 merged:", test_merged_string[2:1:-1])
print("3:1:-1 merged:", test_merged_string[3:1:-1])
print("4:1:-1 merged:", test_merged_string[4:1:-1])
print("5:1:-1 merged:", test_merged_string[5:1:-1])

#
# Strip something from a string
#

stripped = test_merged_string.strip("hel")
print(stripped)

#
# Replace something in a string
#

replaced = test_merged_string.replace("l", "L")
print(replaced)

#
# Using .join() function to a string
#

sample_string = ","
sample_list = ['m','i','a','m','i']
print(sample_string.join(sample_list))