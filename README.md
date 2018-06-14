# SubmitMELA
```
cmsrel CMSSW_9_4_4 
cd CMSSW_9_4_4/src
cmsenv
git clone https://github.com/maravin/SubmitMELA.git
cd SubmitMELA.git
source recipe.sh
scram b -j 8
```

The code really takes only two inputs: inputFile for inputFile and outputFile for output.
For example:
```
produceMELABranches inputFile=DYJets2_11.root outputFile=test.root
```

The code adds the following branches to the tree:
```
Dbkg_VBF, Dbkg_ggH, Dbkg_WH, Dbkg_ZH: discriminators of the VBF, ggH, WH (W->jj), ZH (Z->jj) against Z->tt + 2jets process
ME_*: matrix elements for each signal process
Q2V1, Q2V2, costheta1, costheta2, costhetastar, Phi, Phi1: kinematic inputs to MELA
mjj: invariant mass  of two leading (in pT) jets
```
