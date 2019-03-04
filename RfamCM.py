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

#RNA_final = []
#print(ids[1])
#CM_final = []
#eValue_final = []
#RNA_spare = []
#ids_spare = []
#CM_spare = []
#eValue_spare = []
#spare = []

#for r in RNA:
#	if r in spare:
#		pass
#	else:
#		spare.append(r)
#		n = 0
#		R = []
#		Na = []
#		C = []
#		E = []
#		while n < len(RNA):
#			if r == RNA[n]:
#				R.append(RNA[n])
#				Na.append(ids[n])
#				C.append(CM[n])
#				E.append(eValue[n])
#				n += 1
#			else:
#				n += 1
#		m = 0
#		while m < len(C):
#			value = 0
#			val = 0
#			for v in C:
#				if value == 0:
#					value = v
#					val = v
#					m = m + 1
#					o = 0
#					z = 0
#				else:
#					if float(v) > float(value):
#						value = float(v)
#						o = m
#						m = m + 1
#					elif float(v) < float(val):
#						val = float(v)
#						z = m
#						m = m + 1
#					else:
#						m = m + 1
#
#		RNA_final.append(R[o])
#		ids_final.append(Na[o])
#		CM_final.append(C[o])
#		eValue_final.append(E[o])
#
#		if float(val) != float(value):
#			RNA_spare.append(R[z])
#			ids_spare.append(Na[z])
#			CM_spare.append(C[z])
#			eValue_spare.append(E[z])
#		else:
#			pass

		        #if o == 0:
		#	RNA_spare.append(R[m-1])
		#	names_spare.append(Na[m-1])
		#	CM_spare.append(C[m-1])
		#	eValue_spare.append(E[m-1])
		#else:
		#	RNA_spare.append(R[o-1])
		#	names_spare.append(Na[o-1])
		#	CM_spare.append(C[o-1])
		#	eValue_spare.append(E[o-1])

#print(len(RNA_final))
#print(len(spare))
#print(names_final[1])
#print(len(CM_final))
#print(len(CM_spare))

converter=str(sys.argv[4])

f = open(sys.argv[2], 'r')		#csv file from coordinates2.py
reader = csv.reader(f, delimiter=',')  #1 for Null
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
#print(info[1])

f = open(sys.argv[3], 'w')			#final output
#f.write('Chromosome,Start,End,CM,eValue,Sequence\n')
#k = open(sys.argv[4], 'w')			#spare output
#k.write('Chromosome,Start,End,CM,eValue,Sequence\n')

a = 0
hmm = []
#whyy = []
#b = 0

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
			if name in c: #and RNA_final[a] not in whyy:
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
		i = info[b]
			#nc = i[0]
			#nc = nc.replace(" ",".")
		c = i[0][3:]
		d = i[1]
		e = i[2]
		g = i[3]
			#h = i[5]
			#j = i[6]
		data2 = ''
	        	#if nam in c and g not in hmm:
			#data2 += nc+','
		data2 += c+','
			#data2 += RNA_spare[b]+','
		data2 += d+','
		data2 += e+','
			#data2 += g+','
			#data2 += h+','
		data2 += str(CM[b])+','
		data2 += str(eValue[b])+','
		data2 += g
		f.write(data2)
		f.write('\n')
			        #break
		        #else:
			        #pass
		b = b + 1


f.close()
#print(nc)
#k.close()
print(converter, 'dataset has been generated')
