 
#!/bin/bash

source "$1"

RANDOM_SEED=8095254980

function runREGGAEAssiciationsMode {
	runMethod "REGGAE with association file" reggae "reggae_association_file_new.txt" --scores "${BASE_DIR}/input/deregulated_genes.txt" \
		 --associations "${BASE_DIR}/input/rtis_with_correlation.txt" \
		 --method wrs-test \
		 --bootstrap 0 \
		 --alpha 0.1\
		 --adjust benjamini_yekutieli\
		 --json\
		 --impact pearson_correlation\
		 --confidence-intervals percentile\
		 --sort-rtis-decreasingly \
		 --decreasingly\
		 --abs
}

function runREGGAEMatrixMode {
	runMethod "REGGAE with gene expression matrix" reggae  "reggae_gene_expression_matrix_new.txt" --scores "${BASE_DIR}/input/deregulated_genes.txt" \
		--matrix "${BASE_DIR}/input/matrix.txt" \
		--regulations "${BASE_DIR}/input/rtis.txt" \
		--method wrs-test \
		--bootstrap 0 \
		--alpha 0.1 \
		--adjust benjamini_yekutieli \
		--json \
		--impact pearson_correlation \
		--confidence-intervals percentile \
                --sort-rtis-decreasingly \
		--decreasingly \
		--abs
}

function runMicroREGGAEAssociationMode {
	runMethod "microREGGAE with association file" microReggae  "microReggae_association_mode.txt"\
		--scores "${BASE_DIR}/input/deregulated_genes_microReggae.txt" \
		--associations "${BASE_DIR}/input/rtis_with_correlation_microReggae.txt" \
		--method wrs-test \
		--json 
}


function runMicroREGGAEMatrixMode {
	runMethod "microREGGAE with miRNA expression matrix" microReggae  "microReggae_matrix_mode.txt"\
		--scores "${BASE_DIR}/input/deregulated_genes_microReggae.txt" \
		--matrix "${BASE_DIR}/input/matrix_MicroReggae.txt" \
		--matrix-micro  "${BASE_DIR}/input/microMatrix.txt" \
		--regulations "${BASE_DIR}/input/rtis_microReggae.txt" \
		--method wrs-test \
		--json \
		--impact pearson_correlation 
				
}
function runMicroREGGAEMatrixMode2 {
	runMethod "microREGGAE with miRNA expression matrix2" microReggae  "microReggae_matrix_mode2.txt"\
		--scores "${BASE_DIR}/input/deregulated_genes_microReggae2.txt" \
		--matrix "${BASE_DIR}/input/matrix_MicroReggae.txt" \
		--matrix-micro  "${BASE_DIR}/input/microMatrix.txt" \
		--regulations "${BASE_DIR}/input/rtis_microReggae2.txt" \
		--method wrs-test \
		--json \
		--impact pearson_correlation 
				
}

function runMicroREGGAEMatrixMode3 {
	runMethod "microREGGAE with miRNA expression matrix3" microReggae  "microReggae_matrix_mode3.txt"\
		--scores "${BASE_DIR}/input/deregulated_genes3.txt" \
		--matrix "${BASE_DIR}/input/matrix3.txt" \
		--matrix-micro  "${BASE_DIR}/input/microMatrix3.txt" \
		--regulations "${BASE_DIR}/input/rtis3.txt" \
		--method wrs-test \
		--json \
		--impact pearson_correlation 
				
}
function runMethod {
	mkdir -p "${OUTPUT}"

		
	echo -ne "Running $1 ... "
	"${BINARY_PATH}/$2" --seed "${RANDOM_SEED}" --output "${OUTPUT}/$3" "${@:4}" > /dev/null
	
	if [[ ! -f "${OUTPUT}/$3" ]]; then
		echo "[1;31mNo output has been produced![0m";
		result=1
		return
	fi

	REFERENCE="${BASE_DIR}/results/$3"

	if [[ ! -f ${REFERENCE} ]]; then
		echo "[1;31mCould not open file '${REFERENCE}'[0m";
		result=1
		return
	fi

	local DIFF=$(diff -purN "${REFERENCE}" "${OUTPUT}/$3")

	if [[ "${DIFF}" != '' ]]; then
                echo -e "\e[01;31mError\e[0m"
		echo ""
                echo "    ============================================================"
                echo "    Please type the following commands to check the differences:"
                echo "        > "${BINARY_PATH}/$2" --seed "${RANDOM_SEED}" --output "${OUTPUT}/$3" "${@:4}" > /dev/null" 
                echo "        > vimdiff ${REFERENCE} ${OUTPUT}/$3"
                echo "    ============================================================"
                echo ""
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
#runMicroREGGAEAssociationMode
#runMicroREGGAEMatrixMode
#runMicroREGGAEMatrixMode2
#runMicroREGGAEMatrixMode3
