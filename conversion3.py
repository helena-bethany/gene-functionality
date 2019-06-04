#Convert MAF to stk

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

input_file = sys.argv[1]	#multiz100way multiple alignment file
a = []

alignments = AlignIO.parse(input_file, 'maf')
for alignment in alignments:
	a.append(alignment)	# Extracts each alignment from MAF file

my_records = []
RNA = sys.argv[2]		# Name for output stockholm file
file = open("RNA.stk", 'w')
file.write('# STOCKHOLM 1.0\n')	# Setting up beginning of stk file
file.write('\n')
file.write('#=GF ID '+RNA)
file.write('\n')
file.write('\n')

alignment = a[0]		
for record in alignment:
	n = 1
	ID = record.id		# Extracts ID from alignment information
	seq = str(record.seq)	# Extracts sequence from alignment information
	copies = 0
	while n < len(a):
		other = a[n]
		for record in other:
			d = record.id
			if d == ID:
				seq += str(record.seq)		# MAF file kept the sequence split on numerous lines
				copies += 1			# so this connects them together into one sequence.
			else:
				pass
		n += 1
	ID = ID.replace('.','/')
	if copies != (len(a)-1):
		pass
	else:
		file.write(ID)
		file.write(' '*(40-len(ID)))	#Making sure than spacing is the same for all sequences
		file.write(seq)
		file.write('\n')

file.close()

count = SeqIO.convert("RNA.stk", 'stockholm', RNA+".fa", "fasta")
