#!/bin/bash

source "$1"

RANDOM_SEED=8095254980
ROW_PERMUTATIONS=50000
COL_PERMUTATIONS=5000

function runScoreMethod {
	runMethod "$1" "$2" "$3" --scores "${BASE_DIR}/input/scores.txt"
}

function runIdentifierMethod {
	runMethod "$1" "$2" "$3" --identifier "${BASE_DIR}/input/gse28807.list"
}

function runScoreMethodColwise {
	runMethod "$1" "$2" "$3" \
		--scores "${BASE_DIR}/input/scores_colwise.txt" \
		--pvalue_strategy column-wise \
		--scoring_method independent-shrinkage-t-test \
		--groups "${BASE_DIR}/input/groups.txt" \
		--data_matrix_path "${BASE_DIR}/input/data.txt" \
		--permutations "${COL_PERMUTATIONS}"
}

function runMethod {
	mkdir -p "${OUTPUT}"

	local ID="${1}_${2}"

	echo -n "Running ${ID}... "

	"${BINARY_PATH}/$1" $3 --seed "${RANDOM_SEED}" --categories "${GENERATED_INPUT_DIR}/categories.txt" -o "${OUTPUT}" "${@:4}" > /dev/null

	if [[ ! -f "${OUTPUT}/KEGG.txt" ]]; then
		echo "[1;31mNo output has been produced![0m"
		result=1
		return
	fi

	mv "${OUTPUT}KEGG.txt" "${OUTPUT}${ID}.txt"

	REFERENCE="${BASE_DIR}/results/${ID}.txt"

	if [[ ! -f ${REFERENCE} ]]; then
		echo "[1;31mCould not open file '${REFERENCE}'[0m"
		result=1
		return
	fi

	local DIFF=$(diff -purN "${REFERENCE}" "${OUTPUT}${ID}.txt")

	if [[ "${DIFF}" != '' ]]; then
		echo "[1;31mError[0m"
		result=1
	else
		echo "[0;32mPassed[0m"
		rm "${OUTPUT}${ID}.txt"
		result=0
	fi
}

##
## Run all methods in row-wise mode
##

## enrichment
methods=( sum mean median max-mean )
for i in "${methods[@]}"
do
	runScoreMethod enrichment "${i}_rowwise" "--method $i --permutations ${ROW_PERMUTATIONS}"
	runScoreMethodColwise enrichment "${i}_colwise" "--method $i"
done

runScoreMethod gsea rowwise
runIdentifierMethod gsea identifier_rowwise
runScoreMethodColwise gsea colwise

runScoreMethod weighted-gsea rowwise "--permutations ${ROW_PERMUTATIONS}"
runScoreMethodColwise weighted-gsea colwwise

methods=( 'one-sample-t-test' 'two-sample-t-test' 'two-sample-wilcoxon' )
for i in "${methods[@]}"
do
	runScoreMethod htests "${i}_rowwise" "--method $i"
	runScoreMethodColwise htests "${i}_colwise" "--method $i"
done

runScoreMethod ora rowwise "--reference ${BASE_DIR}/input/reference.txt"
runScoreMethodColwise ora colwise "--reference ${BASE_DIR}/input/reference.txt"
