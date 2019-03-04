#Creating a python script cause it'll be faster than working it out in bash

import sys
import csv

file = open(sys.argv[1], 'r')
lines = file.readlines()

print('CM'+','+'eValue')
for line in lines:
	if '#' in line:
		pass
	else:
		line = " ".join(line.split())
		line = line.replace(" ",',')
		line = line.split(',')
		words = ''
		if line[0] == '1':
			print(line[16]+','+line[17])
		else:
			pass
		#print(words)

file.close()

