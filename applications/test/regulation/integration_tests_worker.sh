#!/bin/bash

source "$1"

RANDOM_SEED=8095254980

function runREGGAEAssiciationsMode {
	runMethod "REGGAE with association file" reggae "reggae_association_file.txt" --scores "${BASE_DIR}/input/deregulated_genes.txt" \
		 --associations "${BASE_DIR}/input/rtis_with_correlation.txt" \
		 --method wrs-test \
		 --bootstrap 0\
		 --alpha 0.1\
		 --adjust benjamini_yekutieli\
		 --json\
		 --impact pearson_correlation\
		 --confidence-intervals percentile\
		 --decreasingly\
		 --abs
}

function runREGGAEMatrixMode {
	runMethod "REGGAE with gene expression matrix" reggae  "reggae_gene_expression_matrix.txt" --scores "${BASE_DIR}/input/deregulated_genes.txt" \
		--matrix "${BASE_DIR}/input/matrix.txt" \
		--regulations "${BASE_DIR}/input/rtis.txt" \
		--method wrs-test \
		--bootstrap 0 \
		--alpha 0.1 \
		--adjust benjamini_yekutieli \
		--json \
		--impact pearson_correlation \
		--confidence-intervals percentile \
		--decreasingly \
		--abs
}

function runMethod {
	mkdir -p "${OUTPUT}"

	echo "Running $1 ..."
	"${BINARY_PATH}/$2" --seed "${RANDOM_SEED}" --output "${OUTPUT}/$3" "${@:4}" > /dev/null
	
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

runREGGAEAssiciationsMode
runREGGAEMatrixMode