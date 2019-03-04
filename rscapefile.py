import sys
import csv

file = open(sys.argv[1], 'r')
x = sys.argv[1]
lines = file.readlines()

for line in lines:
    if "*" in line:
        print(x)
        print(line)
    else:
        pass

file.close()
