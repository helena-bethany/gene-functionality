#Getting chromosome coordinates for RNAcentral dataset

import sys
import csv

IDs = []                                                #IDs is the list of RNAcentral IDs from dataset
names = []						#Names of the ncRNA
sequence = []						#Sequences for each ncRNA
#t = open(sys.argv[3], 'w')				#Reference list of ncRNA from original dataset

with open(sys.argv[1], 'r') as csvFile:			#Converted fasta to csv file from RNAcentral
	reader = csv.reader(csvFile, delimiter = "\t")
	n = 0
	for row in reader:
		sequence.append(row[1])
		info = row[0]
		count = 0
		ID = info[:18]
		IDs.append(ID)
		for letter in info:
			if letter == '(':
				count += 1
			else:
				pass
		if count == 1:
			info = info.split('('  or ')')
			info = info[1]
			info = info.split(')')
			info = info[0]
			info = info.replace(',',';')
		else:
			info = info.split('(')
			num = len(info)
			info = info[num-1]
			info = info.split(')')
			info = info[0]
			info = info.replace(',',';')
		n = n + 1
		names.append(info)
		#t.write(ID)
		#t.write('\n')
		#break

csvFile.close()
print("There were",n,"ncRNA in the dataset.")
#print(IDs[0])
#print(names[0])
#print(sequence[0])

all_coordinates = []                                    #List of containing information from above file

with open('Homo_sapiens.GRCh38.bed') as file:           #HGNC file containing chromosome coordinates for every ncRNA in database
	m = 0
	chrom = []
	chromStart = []
	chromEnd = []
	rcIDs = []
	#databases = []
	for line in file:
		content = []
		content.append(line.strip().split())
		rcID = content[0][3]
		chr = content[0][0]
		chrStart = content[0][1]
		chrEnd = content[0][2]
		#db = content[0][14]
		#db = db.split(',')
		#all = content[0][0:4]
		#all = str(all)
		#all = all.split("'")
		#data = ''
		#allrcID = []
		#for word in all:
		#	if word == '[' or word == ']' or ',' in word:
		#		pass
		#	else:
		#		data += word+','
		if rcID in IDs: #and 'HGNC' in db:
			#all_coordinates.append(data)
			#allrcID.append(rcID)
			chrom.append(chr)
			chromStart.append(chrStart)
			chromEnd.append(chrEnd)
			rcIDs.append(rcID)
			#databases.append(db)
			m = m + 1
			#break
		else:
			pass
			#break

file.close()
print("There were",m,"coordinates found.")
#print(chrom)
#print(chromStart)
#print(chromEnd)
#print(rcIDs)
#print(isinstance(chr,str))
#print(all_coordinates[0])
#print(data)

f = open(sys.argv[2], 'w')                              #File contains location of F and NF ncRNA from dataset							#Count for the dataset info
b = 0							#Count for the coordinates info

f.write('ncRNA,RNAcentral ID,Copies,Chromosome,Start,End,Sequence\n')


for id in rcIDs:
	a = 0						#Count for the dataset info
	for sigh in IDs:
		if id == sigh:
			copies = 0
			for n in rcIDs:
				if sigh == n:
					copies += 1
				else:
					pass
			data = ''
			data += names[a]+','
			data += IDs[a]+','
			data += str(copies)+','
			data += chrom[b]+','
			data += chromStart[b]+','
			data += chromEnd[b]+','
			data += sequence[a]+','
			#data += str(databases[b])
			f.write(data)
			f.write('\n')
			a = a + 1
			#break
		else:
			a = a + 1
			#break
	b = b + 1


f.close()
#print(a)
#print(b)

