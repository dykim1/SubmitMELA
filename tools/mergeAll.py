# code to merge the output from running MELA on svFitted samples by Doyeong. 
import os, re, sys

directory = sys.argv[1]
dir = os.popen("ls " + directory + "/*")
lines = dir.readlines()

datasets = []
for line in lines:
    elements = re.split("-", line)
    file = ""
    if len(re.split("TauTau_13",elements[len(elements)-2])) > 1:
        file = elements[-1]
    else:
        file = elements[-2] + "-" + elements[-1]
    element = re.split("\_", file)[0]
    if not element in datasets:
        print elements
        datasets += [element,]

print datasets


for iSet in datasets:
    #if iSet <> "WJets":
    #    continue

    listOfFiles = []
    dir = os.popen("ls " + directory + "/*" + iSet + "*root")
    files = dir.readlines()
    for iFile in files:
        if len(re.split("-" + iSet + "_", iFile)) > 1:
            listOfFiles += [iFile[:-1],]

    command = "cd " + directory + ";hadd ~/" + iSet + "_svFit_MELA.root " 
    for iFile in listOfFiles:
        newFile = re.split("\/",  iFile)[-1]
        #command  += iFile + " "
        command += newFile + " "
    os.system(command)
    #print command

