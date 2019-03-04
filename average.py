import csv
import sys

total = 0
count = 0
minimum = 0
maximum = 0

switch=str(sys.argv[2])

if switch == "RIsearch":
    with open(sys.argv[1], 'r') as csvFile:
        reader = csv.reader(csvFile, delimiter='\t')
        for line in reader:
            #print(line)
            number = float(line[7])
            total += number
            if number < minimum:
                minimum = number
            elif number > maximum:
                maximum = number
            else:
                pass
            count += 1
    average = float(total/count)
    print(str(minimum)+', '+str(average))
    csvFile.close()

rna=[]
n=0
letters = ['i','E','I']

if switch == "RIblast":
    with open(sys.argv[1], 'r') as File:
        reader = csv.reader(File, delimiter=',')
        for line in reader:
           if line[0][0] in letters:
                pass
           else:
                query=line[1]
                number = float(line[5])
                total += number
                check = n+1
                if len(rna) != check:
                    hits = 1
                    rna.append(query)
                    total += number
                    if number < minimum:
                        minimum = number
                    elif number > maximum:
                        maximum = number
                    else:
                        pass
                    count += 1
                elif query in rna:
                    hits +=1
                    total += number
                    if number < minimum:
                        minimum = number
                    elif number > maximum:
                        maximum = number
                    else:
                        pass
                    count += 1
                else:
                    average= float(total/count)
                    print(rna[n]+', '+str(count)+', '+str(minimum)+', '+str(average))
                    n += 1
                    minimum=0
                    maximum=0
                    total=0
                    count = 0
    File.close()

