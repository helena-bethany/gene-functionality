#Create "random" sequences for the null dataset

import sys
import csv

csv.field_size_limit(sys.maxsize)

start = []
end = []
chromosome = []
length = []

with open(sys.argv[1], 'r') as csvFile:			#final csv file with sequences
	reader = csv.reader(csvFile, delimiter = ",")
	for row in reader:
		s=0
		e=0
		if 'C' in row[0]:
			pass
		else:
			chromosome.append(row[0])
		if 'S' in row[1]:
			pass
		else:
			s += int(row[1])
			start.append(s)
		if 'E' in row[2]:
			pass
		else:
			e += int(row[2])
			end.append(e)
		l = e-s
		length.append(l)


csvFile.close()
length = length[1:]
start = start[1:]
end = end[1:]
chromosome = chromosome[1:]

s = open(sys.argv[2], 'w')
s.write('Chromosome,Start,End,Sequence\n')
e = open(sys.argv[3], 'w')
e.write('Chromosome,Start,End,Sequence\n')
n = 1
total = 0

with open('GRCh38_genome.csv', 'r') as csvFile:                 #Converted fasta to csv file from chromosome
	reader = csv.reader(csvFile, delimiter = "\t")
	chr = []
	for row in reader:
		if 'NC' in row[0]:
			#print(row[0])
			chr.append(str(row[1]))
			total += 1
		else:
			pass
		
		

while n <= total:
	a = 0
	for chromo in chromosome:
		if chromo is 'X':
			chromo = 23
		elif chromo is 'Y':
			chromo = 24
		elif chromo is 'M':
			chromo = 25
		else:
			pass
		if int(chromo) == n:
			z = n - 1
			sequence_start = ''
			right1 = int(start[a]) - 20001
			left1 = right1 - int(length[a])
			sequence_start += ('Chr'+str(chromosome[a]))+','
			sequence_start += str(left1)+','
			sequence_start += str(right1)+','
			seq = chr[z][left1:right1]
			sequence_start += seq
			s.write(sequence_start)
			s.write('\n')
			sequence_end = ''
			left2 = int(end[a]) + 19999
			right2 = left2 + int(length[a])
			sequence_end += ('Chr'+str(chromosome[a]))+','
			sequence_end += str(left2)+','
			sequence_end += str(right2)+','
			seq = chr[z][left2:right2]
			sequence_end += seq
			e.write(sequence_end)
			e.write('\n')
		else:
			pass
		a = a + 1
	n = n + 1



csvFile.close()
print('Null dataset sequences generated')
