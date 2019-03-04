import sys
import csv

info1 = []
info2 = []

with open("intact.txt", 'r') as csvFile:
    reader = csv.reader(csvFile, delimiter = '\t')
    for line in reader:
        if "rnacentral" in line[2] or "rnacentral" in line[3]:
            string = ''
            string += line[0]+','
            string += line[1]+','
            string += line[2]+','
            string += line[3]+','
            string += line[11]+','
            string += line[14]+','
            info1.append(string)
        else:
            pass

      
with open("Homo_sapiens.GRCh38.bed", 'r') as File:
    reader = csv.reader(File, delimiter = '\t')
    for line in reader:
        string = ''
        string += line[0]+','
        string += line[1]+','
        string += line[2]+','
        string += line[3]
        info2.append(string)
    
csvFile.close()
File.close()

print("Total RNAcentral in INTACT is",len(info1))
print("Total chormosome coordinates is",len(info2))


output = open("181206intactlocations.csv",'w')
output.write("rnacentral,chromosome,start,end,interactorA,interactorB,alternative_interactorA,alternative_interactorB,interactionType,confidenceValues")
x = 0
count = 0

while x < len(info2):
    location = info2[x]
    location = location.split(',')
    chromo = location[0]
    start = location[1]
    end = location[2]
    rnacent = location[3]
    #print(location)
    y = 0
    while y < len(info1):
        rna = info1[y]
        rna = rna.split(',')
        #print(rna)
        if rnacent in rna[2] or rnacent in rna[3]:
            string = ''
            string += rnacent+','
            string += chromo+','
            string += start+','
            string += end+','
            string += rna[0]+','
            string += rna[1]+','
            string += rna[2]+','
            string += rna[3]+','
            string += rna[4]+','
            string += rna[5]
            output.write(string)
            output.write('\n')
            y += 1
            count += 1
        else:
            y += 1
    x += 1

output.close()
print(count,"coordinates converted")
