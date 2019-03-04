import sys
import csv

csv.field_size_limit(sys.maxsize)

chromo = []
start = []
end = []

with open(sys.argv[1], 'r') as File:
    reader = csv.reader(File, delimiter = ',')
    for line in reader:
        if 'Start' in line[2]:
            pass
        else:
            chromo.append(str(line[1]))
            start.append(str(line[2]))
            end.append(str(line[3]))

File.close()
#print(chromo[1])
#print(start[1])
#print(end[1])
#print(len(chromo))

file = open(sys.argv[3], 'w')
file.write("Chromosome,Start,End,Expression,Average\n")
n = 0

with open(sys.argv[2], 'r') as csvFile:
    reader = csv.reader(csvFile, delimiter = '\t')
    while n < len(chromo):
        c = 'chr'+chromo[n]
        s = int(start[n])
        e = int(end[n])
        exp = 0
        for line in reader:
            if c == line[0]:
                #print(line[1])
                #print(line[2])
                if int(line[1]) == s:
                    exp += float(line[3])
                elif int(line[1]) < s < int(line[2]):
                    exp += float(line[3])
                elif s < int(line[1]) < e and s < int(line[2]) < e:
                    exp += float(line[3])
                elif int(line[1]) < e < int(line[2]):
                    exp += float(line[3])
                elif int(line[2]) == e:
                    exp += float(line[3])
                elif int(line[1]) < s:
                    pass
                elif int(line[2]) > e:
                    break
                else:
                    pass
            else:
                pass
        seq = e-s
        average = float(exp/seq)
        output=c+','+str(s)+','+str(e)+','+str(exp)+','+str(average)
        file.write(output)
        file.write('\n')
        n += 1
        
csvFile.close()
