pushd $CMSSW_BASE/src

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement ; git checkout -b v2.1.7b1 ; . setup.sh -j 12)
cp "$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/"*so "$CMSSW_BASE/lib/$SCRAM_ARCH/"
popd