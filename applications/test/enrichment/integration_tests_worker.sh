source $1

RANDOM_SEED=8095254980
PERMUTATIONS=50000

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
		--groups ${BASE_DIR}/input/groups.txt \
		--data_matrix_path ${BASE_DIR}/input/data.txt \
		--permutations ${PERMUTATIONS}
}

function runMethod {
	mkdir -p "${OUTPUT}"

	local ID="${1}_${2}"

	echo -n "Running ${ID}... "

	"${BINARY_PATH}/$1" $3 --seed ${RANDOM_SEED} --categories "${GENERATED_INPUT_DIR}/categories.txt" -o "${OUTPUT}" ${@:4} > /dev/null

	if [[ ! -f "${OUTPUT}/KEGG.txt" ]]; then
		echo "No output has been produced!"
		exit 1
	fi

	mv "${OUTPUT}KEGG.txt" "${OUTPUT}${ID}.txt"

	REFERENCE="${BASE_DIR}/results/${ID}.txt"

	if [[ ! -f ${REFERENCE} ]]; then
		echo "Could not open file '${REFERENCE}'"
		exit 1
	fi

	local DIFF=`diff -purN "${REFERENCE}" "${OUTPUT}${ID}.txt"`

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
for i in ${methods[@]}
do
	runScoreMethod enrichment "${i}_rowwise" "--method $i --permutations ${PERMUTATIONS}"
done

runScoreMethod gsea rowwise
runIdentifierMethod gsea identifier_rowwise
runScoreMethodColwise gsea colwise
runIdentifierMethodColwise gsea identifier_colwise

runScoreMethod weighted-gsea rowwise "--permutations ${PERMUTATIONS}"

methods=( 'one-sample-t-test' 'two-sample-t-test' 'two-sample-wilcoxon' )
for i in ${methods[@]}
do
	runScoreMethod htests "${i}_rowwise" "--method $i"
done

runScoreMethod ora rowwise "--reference ${BASE_DIR}/input/reference.txt"

