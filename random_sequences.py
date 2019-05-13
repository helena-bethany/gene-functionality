#Create "random" sequences for the null dataset

import sys
import csv

csv.field_size_limit(sys.maxsize)

start = []
end = []
#chromosome = []

with open(sys.argv[1], 'r') as csvFile:
	reader = csv.reader(csvFile, delimiter = ',')
	for row in reader:
		#chromosome.append(row[3])
		if 'S' in row[4]:
			pass
		else:
			no = int(row[4])
			start.append(no)
		if 'E' in row[5]:
			pass
		else:
			mo = int(row[5])
			end.append(mo)


csvFile.close()
#print(chromosome)
#print(start[0])
#print(end)

chr = ''
s = open('random_left.csv', 'w')
e = open('random_right.csv', 'w')

with open(sys.argv[2], 'r') as csvFile:                 #Converted fasta to csv file from chromosome
	reader = csv.reader(csvFile, delimiter = "\t")
	for row in reader:
		chr += str(row[1])
		break
	count = 0
	while count <= len(chr):
		for number in end:
			if count == number:
				sequence = ''
				left = count + 20000
				right = left + 150
				sequence += chr[left:right]
				e.write(sequence)
				e.write('\n')
				break
			else:
				pass
		for number in start:
			if count == number:
				sequence = ''
				right = count - 20000
				left = right - 150		#150 is the average snRNA length
				sequence += chr[left:right]
				s.write(sequence)
				s.write('\n')
				break
			else:
				pass

		count = count + 1
		#counter = count/1000000
		#if isinstance(counter,int) == 'TRUE':
		#	print(count)
		#else:
		#	pass

		#count = count + 1

csvFile.close()
#print(len(chr))
#print(count)
