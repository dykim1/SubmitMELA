import os, sys, re

dirToCheck = os.popen("ls -lrt " + sys.argv[1])
files = dirToCheck.readlines()

for file in files:
    elements = re.split(" ", file)
    # check if it is a header or total or something like that
    if len(elements) < 4:
        continue
    counter = 4
    toCheck = elements[counter]
    while elements[counter] == '':
        counter +=  1
        toCheck = elements[counter]
        
    mem = int(toCheck)
    if mem < 10000:
        print file
