source $1

function runMethod {
	mkdir -p "${OUTPUT}"

	local ID="${1}_${2}"

	echo -n "Running ${ID}... "

	"${BINARY_PATH}/$1" $3 --seed ${RANDOM_SEED} --scores "${BASE_DIR}/input/scores.txt" --categories "${GENERATED_INPUT_DIR}/categories.txt" -o "${OUTPUT}" > /dev/null

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
	runMethod enrichment "${i}_rowwise" "--method $i --permutations ${PERMUTATIONS}"
done

runMethod gsea rowwise
runMethod weighted-gsea rowwise "--permutations ${PERMUTATIONS}"

methods=( 'one-sample-t-test' 'two-sample-t-test' 'two-sample-wilcoxon' )
for i in ${methods[@]}
do
	runMethod htests "${i}_rowwise" "--method $i"
done

runMethod ora rowwise "--reference ${BASE_DIR}/input/reference.txt"

