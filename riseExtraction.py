import sys
import csv

#file1 = open(sys.argv[1], 'r')        #rise dataset
#lines = file1.readlines()

coordinates = []

with open(sys.argv[2], 'r') as csvFile:
    reader = csv.reader(csvFile, delimiter = ',')
    for line in reader:
        if "Start" in line:
            pass
        else:
            chromo = line[4]
            start = line[5]
            end = line[6]
            coords = "chr"+str(chromo)+":"+str(start)+":"+str(end)
            coordinates.append(coords)

csvFile.close()

output = open(sys.argv[3], 'w')
total = 0
firstRNA = 0
secondRNA = 0
found = 0

with open(sys.argv[1], 'r') as file:
    reader = csv.reader(file, delimiter = '\t')
    for line in reader:
        if "." in line[0] or "." in line[1] or "." in line[2]:
            chromo1 = "chr0"
            start1 = 0
            end1 = 0
        else:
            chromo1 = line[0]
            #print(chromo1)
            start1 = int(line[1])
            end1 = int(line[2])
            name1 = line[14]
            firstRNA += 1
            #coords1 = chromo1+":"+str(start1)+"-"+str(end1)
            #print(coords1)
        if "." in line[3] or "." in line[4] or "." in line[5]:
            chromo2 = "chr0"
            start2 = 0
            end2 = 0
        else:
            chromo2 = line[3]
            start2 = int(line[4])
            end2 = int(line[5])
            #coords2 = chromo2+":"+str(start2)+"-"+str(end2)
            name2 = line[15]
            secondRNA += 1
            #print(coords2)
        for rna in coordinates:
            rn = rna.split(":")
            c = rn[0]
            #print(c)
            s = int(rn[1])
            e = int(rn[2])
            count = 0
            interact = []
            if chromo1 in c:
                if start1 <= s <= end1 or start1 <= e <= end1:
                    interact.append(name1)
                    count += 1
                elif start2 <= s <= end2 or start2 <= e <= end2:
                    interact.append(name2)
                    count += 1
                else:
                    pass
            else:
                pass
            if count > 0:
                string = ''
                string += rna+','+str(count)+','
                for x in interact:
                    string += x+','
                output.write(string)
                output.write("\n")
                found += 1
            else:
                pass
        total += 1
        #break

file.close()
print("Total of RNA1 is",firstRNA) 
print("Total of RNA2 is",secondRNA)
print("Total number of RNA is",total)
print("Matches found is",found)

#print(coordinates[1])
