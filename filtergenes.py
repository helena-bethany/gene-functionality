import sys
import csv

csv.field_size_limit(sys.maxsize)

locations = []

with open(sys.argv[1], 'r') as File:
    reader = csv.reader(File, delimiter = '\t')
    for line in reader:
        n = 0
        while n < len(line):
            content = line[n]
            if 'NC_' in content:
                string = ''
                string+=str(line[n-1])+":"+str(line[n+1])+":"+str(line[n+2])
                locations.append(string)
            else:
                pass
            n += 1
            
File.close()
location = locations[:250]

output = open(sys.argv[2], "w")
output.write("Chr,Start,End,Sequence\n")

with open("/home/helenacooper/ncRNA_Project/GRCh38_genome.csv", 'r') as csvFile:
    reader = csv.reader(csvFile, delimiter = '\t')
    chr = []
    for row in reader:
        if 'NC' in row[0]:
            chr.append(str(row[1]))
        else:
            pass

csvFile.close()

for where in location:
    where = where.split(":")
    chromo = where[0]
    start = int(where[1])
    end = int(where[2])
    n = 1
    if chromo is 'X':
        chrom = 23
    elif chromo is 'Y':
        chrom = 24
    elif chromo in 'MT':
        chrom = 25
    else:
        chrom = int(chromo)
    while n <= len(chr):
        if chrom == n:
            z = n - 1
            seq = chr[z][start:end]
            sequence = ''
            sequence += chromo+','+str(start)+','+str(end)+','+seq
            output.write(sequence)
            output.write('\n')
            n = len(chr)
        else:
            n += 1

output.close()
