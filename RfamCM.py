#Create files with the largest CM score

import sys
import csv

f = open(sys.argv[1], 'r')			#text file from Rfam batch search
lines = f.readlines()

RNA = []
ids = []
CM = []
eValue = []

for line in lines:
	if '#' in line:
		pass
	else:
		line = " ".join(line.split())
		line = line.replace(" ",',')
		line = line.split(',')
		#words = ''
		if line[0] == '1':
			RNA.append(line[1])
			ids.append(line[3])
			CM.append(line[16])
			eValue.append(line[17])
		else:
			pass

converter=str(sys.argv[4])

f = open(sys.argv[2], 'r')		#csv file from coordinates2.py
reader = csv.reader(f, delimiter=',')   #1 for Null
copy = []
info = []

if converter == 'Test':
	for row in reader:
		if row[3] == "Start":
			pass
		else:
			info.append(row)
if converter == 'Null':
	for row in reader:
		if row[1] == "Start":
			pass
		else:
			info.append(row)

f.close()

f = open(sys.argv[3], 'w')			#final output
#f.write('Chromosome,Start,End,CM,eValue,Sequence\n')
#k = open(sys.argv[4], 'w')			#spare output
#k.write('Chromosome,Start,End,CM,eValue,Sequence\n')

a = 0
hmm = []

if converter == 'Test':
	f.write('Chromosome,Start,End,CM,eValue,Sequence\n')
	while a < len(ids):
		name = ids[a]
	        #nam = ids_spare[a]
		for i in info:
			nc = i[0]
		        #nc = nc.replace(" ",".")
			c = i[1]
			d = i[2]	#chromosome
			e = i[3]	#start
			g = i[4]	#end
			h = i[5]	#sequence
			data = ''
			if name in c:
			        #data += nc+','
			        #data += c+','
			        #data += RNA_final[a]+','
				data += d+','
				data += e+','
				data += g+','
				data += str(CM[a])+','
				data += str(eValue[a])+','
				data += h
				f.write(data)
				f.write('\n')
				hmm.append(g)
			        #why.append(RNA_final[a])
				break
			else:
				pass
		a = a + 1

b = 0

if converter == 'Null':
	while b < len(ids):
		name = ids[b]
		for i in info:
			c = i[0][3:]
			d = i[1]
			e = i[2]
			g = i[3]
			if name in d:
				data2 = ''
				data2 += c+','
				data2 += d+','
				data2 += e+','
				data2 += str(CM[b])+','
				data2 += str(eValue[b])+','
				data2 += g
				f.write(data2)
				f.write('\n')
				break
			else:
				pass
		b = b + 1


f.close()
print(converter, 'dataset has been generated')

