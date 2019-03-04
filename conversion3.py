#biopython to the rescue

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

input_file = sys.argv[1]
a = []

alignments = AlignIO.parse(input_file, 'maf')
for alignment in alignments:
	a.append(alignment)

#print(len(a))

#output_file = 'maf'.format(os.path.splitext(input_file)[0])
#f = open("output_file.maf", 'w')
my_records = []
RNA = sys.argv[2]
file = open("RNA.stk", 'w')
file.write('# STOCKHOLM 1.0\n')
file.write('\n')
file.write('#=GF ID '+RNA)
file.write('\n')
file.write('\n')

alignment = a[0]
for record in alignment:
	n = 1
	ID = record.id
	seq = str(record.seq)
	copies = 0
	while n < len(a):
		other = a[n]
		for record in other:
			d = record.id
			if d == ID:
				seq += str(record.seq)
				copies += 1
			else:
				pass
		n += 1
	ID = ID.replace('.','/')
	#rec = SeqRecord(Seq(seq), id=ID)
	#my_records.append(rec)
	#my_records = my_records.replace('.','/')
	if copies != (len(a)-1):
		pass
	else:
		file.write(ID)
		file.write(' '*(40-len(ID)))
		file.write(seq)
		file.write('\n')
#SeqIO.write(my_records, "output_file.maf", 'maf')


#print(seq)

file.close()

count = SeqIO.convert("RNA.stk", 'stockholm', RNA+".fa", "fasta")
#print("Converted %i records" %count)


#records = SeqIO.parse(input_file, 'maf')
#records = list(records)

#maxlen = max(len(record.seq) for record in my_records)
#print(maxlen)

#for record in my_records:
#	if len(record.seq) != maxlen:
#		sequence = str(record.seq).ljust(maxlen, '.')
#		record.seq = Seq(sequence)
#assert all(len(record.seq) == maxlen for record in my_records)

#output_file = 'maf'.format(os.path.splitext(input_file)[0])
#with open("RNA.stk", 'w') as f:
#	SeqIO.write(my_records, f, 'stockholm')
#alignment = AlignIO.read(output_file, 'maf')
#print(alignment)

#alignments = AlignIO.parse(output_file, 'maf')
#count = AlignIO.write(alignments, "RNA.stk", "stockholm")
#print("Converted %i records" %count)

#idx = AlignIO.MafIO.MafIndex("mafOut.mafindex", "mafOut", "hg38")



