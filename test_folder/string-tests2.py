#
# Code in this file is based on a Medium article: https://towardsdatascience.com/working-with-strings-in-python-a-cheat-sheet-e6f462e611f0
# Implemented and written by Zoltan Dul, DMD, Phd
#

#
# Using Pandas
#

import pandas as pd
s = pd.Series(["Hello, world"])

#
# Counting Characters
#

sent = 'His name is Abdul. Abdul is a good boy. Boy o boy'
counted_chars = str(sent.count('a'))
print("I counted " + counted_chars + " chars of \"a\".")

#
# Counting the occurance of an word
#

counted_chars2 = str(sent.count('name'))
print("I counted " + counted_chars + " chars of \"name\".")

#
# Count ALL characters in a string using Counter from Collections lib.
#

from collections import Counter

counting_the_chars = Counter('aaaaabbbbbbxxxxxxdzzzzffff')
print(counting_the_chars)

#
# Counting ALL words in a sentence based on the previous lib
#

sent = 'His name is Abdul. Abdul is a good boy, boy o boy'
word_list = sent.split(' ')

# Count all worlds
print(Counter(word_list))

# Count the most common 3
print(Counter(word_list).most_common(3))

#
# Access the index of a character/word
#

sample_string = "I like to speak loud, and enjoy my life"
print("Position of word \"loud\" stars at " + str(sample_string.find("loud")))
print("Character at 5th position is \"" + sample_string[5] + "\".")

#
# Using F-Strings
# value:width.precision
#

name = 'Zoltan'
age = 33
print(f"my name is {name} and I am {age} years old")

#
# Using RegEx
#

import re

text = 'hello my phone number is (000)-111-2222'
pattern = 'phone'
ret1 = re.search(pattern, text)
print(ret1)

text = 'my zipcode is 22202'
pattern = r'\d\d\d\d\d'
ret2 = re.search(pattern, text)
print(ret2)

text = 'the dog and the hog sat in the fog'
ret3 = re.findall(r'.og', text)
print(ret3)
