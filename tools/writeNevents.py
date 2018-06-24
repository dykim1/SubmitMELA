import os, sys, re

#For every root file in origin directory, gets nevents histogram; gets tree from a similar file in destination, and writes tree/nevents combo in a local directory
#Example 
#    python writeNevents.py /scratch/kkaadze/SkimmedNtuples_SMHTT2016/smhmt_20march/ /hdfs/store/user/ymaravin/2016/mt/ mutau_tree

if len(sys.argv) != 4:
    print "Need 3 arguments: python writeNevents.py originDir destinationDir treeName"
    sys.exit(0)



origin = sys.argv[1]
destination =  sys.argv[2]
treeName = sys.argv[3]

dir = os.popen("ls " + origin + "/*root")
ifiles = dir.readlines()

for i in ifiles:
    fileName = re.split("\/", i)[-1][:-1]
    #check that an origin file exists in destination
    dir = os.popen("ls "  + destination + "/" + fileName[:-5] + "_svFit_MELA.root")
    if len(dir.readlines()) == 0:
        continue
    print fileName
    os.system("cp " + destination + "/" + fileName[:-5] + "_svFit_MELA.root .")
    file = open("macro.C", "w")
    file.write("void macro() {\n")
    file.write("    TFile* in2 = TFile::Open(\"" + i[:-1] + "\");\n")
    file.write("    TH1F* h = (TH1F*)in2->Get(\"nevents\");\n")
    file.write("    TFile* in1 = TFile::Open(\"" + fileName[:-5] + "_svFit_MELA.root\", \"update\");\n")
    file.write("    h->Write();\n")
    file.write("    in2->Close();\n")
    file.write("    in1->Close();\n")
    file.write("    gROOT->ProcessLine(\".q\");\n")
    file.write("}\n")
    file.close()
    os.system("root -l -b macro.C")
    os.system("rm macro.C")


