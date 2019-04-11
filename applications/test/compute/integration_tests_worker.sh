#!/bin/bash

source "$1"


function runScoreMode {
	runMethod "computeNetworkFormular with association file" computeNetworkFormular "resultsFoldChange.txt" --scores "${BASE_DIR}/input/geneScores.txt" \
		 --method log-mean-fold-quotient \
		 --scores-mirna "${BASE_DIR}/input/miRNAScoresFoldChange.txt" \
		 --regulations "${BASE_DIR}/input/targets.txt" 
}

function runScoreModeTTest {
	runMethod "computeNetworkFormular with association file - independent t-test" computeNetworkFormular "resultsIndependentTTest.txt" --scores "${BASE_DIR}/input/geneScores2.txt" \
		 --method independent-t-test \
		 --scores-mirna "${BASE_DIR}/input/miRNAScores_independent-t-test-belongsTo2.txt" \
		 --regulations "${BASE_DIR}/input/targets2.txt" 
}
function runScoreModeShrinkage {
	runMethod "computeNetworkFormular with association file - independent shrinkage t-test" computeNetworkFormular "resultsIndependentShrinkageTTest.txt" --scores "${BASE_DIR}/input/geneScores2.txt" \
		 --method independent-shrinkage-t-test \
		 --scores-mirna "${BASE_DIR}/input/miRNAScores_shrinkage-t-test-belongsTo2.txt" \
		 --regulations "${BASE_DIR}/input/targets2.txt" 
}

function runMatrixMode {
	runMethod "computeNetworkFormular with gene expression matrix- foldChange" computeNetworkFormular  "resultsFoldChange.txt" 	--scores "${BASE_DIR}/input/geneScores.txt" \
		--expression-matrix "${BASE_DIR}/input/miRNAmatrix.txt" \
		--regulations "${BASE_DIR}/input/targets.txt" \
		--method log-mean-fold-quotient \
		--groups "${BASE_DIR}/input/groups.txt" 
		
}
function runMatrixModeTTest {
	runMethod "computeNetworkFormular with gene expression matrix- independent t-test" computeNetworkFormular  "resultsIndependentTTest.txt" 	--scores "${BASE_DIR}/input/geneScores2.txt" \
		--expression-matrix "${BASE_DIR}/input/miRNAmatrix2.txt" \
		--regulations "${BASE_DIR}/input/targets2.txt" \
		--method independent-t-test \
		--groups "${BASE_DIR}/input/groups2.txt" 
		
}

function runMatrixModeShrinkageTTest {
	runMethod "computeNetworkFormular with gene expression matrix- independent shrinkage t-test" computeNetworkFormular  "resultsIndependentShrinkageTTest.txt" 	--scores "${BASE_DIR}/input/geneScores2.txt" \
		--expression-matrix "${BASE_DIR}/input/miRNAmatrix2.txt" \
		--regulations "${BASE_DIR}/input/targets2.txt" \
		--method independent-shrinkage-t-test \
		--groups "${BASE_DIR}/input/groups2.txt" 
		
}

function runNOD {
    runMethod "NOD" computeNOD "resultsNOD.txt" --scores "${BASE_DIR}/input/geneScores_NOD.txt" \
    --matrix  "${BASE_DIR}/input/mRNAMatrix.txt"\
    --matrix-micro "${BASE_DIR}/input/miRNAmatrix_NOD.txt"\
    --regulations "${BASE_DIR}/input/targets_NOD.txt"\
    --scoresMirna "${BASE_DIR}/input/scoresMirna_NOD.txt"\
    --treshold -0.0
}


function runNOD2 {
    runMethod "NOD2" computeNOD "resultsNOD2.txt" --scores "${BASE_DIR}/input/geneScores_NOD2.txt" \
    --matrix  "${BASE_DIR}/input/mRNAMatrix.txt"\
    --matrix-micro "${BASE_DIR}/input/miRNAmatrix_NOD.txt"\
    --regulations "${BASE_DIR}/input/targets_NOD.txt"\
    --scoresMirna "${BASE_DIR}/input/scoresMirna2_NOD.txt"\
    --treshold -0.3
}

function runNOD3 {
    runMethod "NOD3" computeNOD "resultsNOD3.txt" --scores "${BASE_DIR}/input/geneScores_NOD2.txt" \
    --matrix  "${BASE_DIR}/input/mRNAMatrix_NOD3.txt"\
    --matrix-micro "${BASE_DIR}/input/miRNAmatrix_NOD.txt"\
    --regulations "${BASE_DIR}/input/targets_NOD.txt"\
    --scoresMirna "${BASE_DIR}/input/scoresMirna3_NOD.txt"\
    --treshold -0.6
}

function runMethod {
	mkdir -p "${OUTPUT}"
	echo "Running $1 ..."
	"${BINARY_PATH}/$2" --output "${OUTPUT}/$3" "${@:4}" > /dev/null
	
	if [[ ! -f "${OUTPUT}/$3" ]]; then
		echo "[1;31mNo output has been produced![0m"
		result=1
		return
	fi

	REFERENCE="${BASE_DIR}/results/$3"

	if [[ ! -f ${REFERENCE} ]]; then
		echo "[1;31mCould not open file '${REFERENCE}'[0m"
		result=1
		return
	fi

	local DIFF=$(diff -purN "${REFERENCE}" "${OUTPUT}/$3")

	if [[ "${DIFF}" != '' ]]; then
		echo "[1;31mError[0m - Type to check: vimdiff ${REFERENCE} ${OUTPUT}/$3"
		result=1
	else
		echo "[0;32mPassed[0m"
		rm "${OUTPUT}/$3"
		result=0
	fi
}

##
## Run all methods
##
runScoreMode
runMatrixMode
runMatrixModeTTest
runMatrixModeShrinkageTTest
runScoreModeTTest
runScoreModeShrinkage
runNOD
runNOD2
runNOD3