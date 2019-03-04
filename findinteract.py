import sys
import csv

for open("181206intactlocations.csv", 'r') as csvFile:
    reader = csv.reader(csvFile, delimiter = ',')
