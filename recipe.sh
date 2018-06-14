pushd $CMSSW_BASE/src

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement ; git checkout -b v2.1.7b1 ; . setup.sh -j 12)
popd