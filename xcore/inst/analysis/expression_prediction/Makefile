PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
NCPU = 6
IMAGE_BIN = $(PROJECT_DIR)/IMAGE/IMAGE.pl
TEPIC_BIN = $(PROJECT_DIR)/TEPIC/Code/TEPIC.sh
INVOKE_DIR = $(PROJECT_DIR)/TEPIC/MachineLearningPipelines/INVOKE/
INVOKE_OUT = TEPIC_out/INVOKE/

.PHONY: clean
clean: ## delete directories created by runing 'image', 'tepic_invoke'
	@echo "cleaning up"
	rm -r IMAGE_out ${INVOKE_OUT}

.PHONY: install
install: ## download programs and set up environment
	@echo "setting up IMAGE"
	wget -c http://bioinformatik.sdu.dk/solexa/webshare/IMAGE/IMAGE_v1.1.tar.gz
	tar -xzf IMAGE_v1.1.tar.gz
	wget -c -P IMAGE/ https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
	chmod a+x IMAGE/IMAGE.pl
	gunzip IMAGE/mm9.fa.gz
	@echo "setting up TEPIC"
	git clone https://github.com/SchulzLab/TEPIC.git
	cd ./TEPIC/Code && ./compile_TRAP_install_R_packages.sh && cd ../..
	wget -c -P ${INVOKE_DIR}/ExampleData/ http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
	gunzip ${INVOKE_DIR}/ExampleData/chr1.fa.gz
	wget -c -P ${INVOKE_DIR}/ExampleData/ http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
	gunzip ${INVOKE_DIR}/ExampleData/gencode.v26.annotation.gtf.gz

.PHONY: image
image: ## run IMAGE on it's example data
	@echo "Running IMAGE"
	mkdir IMAGE_out
	${IMAGE_BIN} \
		-region IMAGE/examples/Enhancers.txt \
		-expression IMAGE/examples/GeneExpression.txt \
		-fasta IMAGE/mm9.fa \
		-RNADesign 1 1 2 2 3 3 4 4 \
		-EnhancerDesign 1 1 2 2 3 3 4 4 \
		-p ${NCPU} \
		-n IMAGE_out/IMAGE_test
	mv IMAGE_test.R IMAGE_test.Rdata

.PHONY: tepic_invoke
tepic_invoke: ## run TEPIC INVOKE pipeline
	@echo "Running TEPIC INVOKE"
	mkdir -p ${INVOKE_OUT}
	${TEPIC_BIN} \
		-g ${INVOKE_DIR}/ExampleData/Homo_sapiens.GRCh38.dna.chromosome.1.fa \
		-b ${INVOKE_DIR}/ExampleData/S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1.bed \
		-o ${INVOKE_OUT}/INVOKE_test \
		-p TEPIC/PWMs/1.0/pwm_vertebrates_jaspar_uniprobe_original.PSEM \
		-a ${INVOKE_DIR}/ExampleData/gencode.v26.annotation.gtf.gz \
		-w 3000 \
		-e FALSE
	python ${INVOKE_DIR}/Scripts/integrateData.py ${INVOKE_OUT}/*Peak_Features_Affinity_Gene_View_Filtered.txt ${INVOKE_DIR}/ExampleData/S001S712_gene_quantification_rsem_grape2_crg_GRCh38_20150622_chr1.txt ${INVOKE_OUT}/IntegratedData/input.txt
	Rscript Scripts/INVOKE.R \
		--dataDir=${INVOKE_OUT}/IntegratedData/ \
		--outDir=${INVOKE_OUT} --response=Expression \
		--cores ${NCPU} \
		--alphas=0.1 \
		--regularisation=E \
		--testsize=0.2 \
		--innerCV=6 \
		--outerCV=3 \
		--performance=TRUE \
		--ftest=FALSE


# Based on http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.DEFAULT_GOAL := help
.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
