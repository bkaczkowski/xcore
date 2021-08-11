PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
NCPU = 6
IMAGE_DIR = $(PROJECT_DIR)/IMAGE
IMAGE_BIN = $(IMAGE_DIR)/IMAGE.pl
IMAGE_OUT = IMAGE_out
TEPIC_DIR = $(PROJECT_DIR)/TEPIC
TEPIC_BIN = $(TEPIC_DIR)/Code/TEPIC.sh
TEPIC_OUT = TEPIC_out
INVOKE_DIR = $(TEPIC_DIR)/MachineLearningPipelines/INVOKE
INVOKE_BIN = $(INVOKE_DIR)/runInvokeAnalysis.sh
INVOKE_OUT = ${TEPIC_OUT}/INVOKE
DYNAMITE_DIR = $(TEPIC_DIR)/MachineLearningPipelines/DYNAMITE
DYNAMITE_BIN = $(DYNAMITE_DIR)/runDYNAMITE.sh
DYNAMITE_OUT = ${TEPIC_OUT}/DYNAMITE

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
	chmod a+x ${TEPIC_BIN} ${INVOKE_BIN} ${DYNAMITE_BIN}
	wget -c -P ${INVOKE_DIR}/ExampleData/ http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
	gunzip ${INVOKE_DIR}/ExampleData/chr1.fa.gz
	wget -c -P ${INVOKE_DIR}/ExampleData/ http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
	gunzip ${INVOKE_DIR}/ExampleData/gencode.v26.annotation.gtf.gz

.PHONY: image
image: ## run IMAGE on it's example data
	@echo "Running IMAGE"
	mkdir -p ${IMAGE_OUT}
	${IMAGE_BIN} \
		-region ${IMAGE_DIR}/examples/Enhancers.txt \
		-expression ${IMAGE_DIR}/examples/GeneExpression.txt \
		-fasta ${IMAGE_DIR}/mm9.fa \
		-RNADesign 1 1 2 2 3 3 4 4 \
		-EnhancerDesign 1 1 2 2 3 3 4 4 \
		-p ${NCPU} \
		-n ${IMAGE_OUT}/IMAGE_test
	mv ${IMAGE_OUT}/IMAGE_test.R ${IMAGE_OUT}/IMAGE_test.Rdata

# TRAP takes 41s on 105K BED file with 191K genome
.PHONY: tepic
tepic: ## run TEPIC on it's example data
	@echo "Running TEPIC INVOKE"
	mkdir -p ${TEPIC_OUT}
	${TEPIC_BIN} \
		-g ${TEPIC_DIR}/Test/example_sequence.fa \
 		-b ${TEPIC_DIR}/Test/example_regions.bed \
 		-o ${TEPIC_OUT}/TEPIC-Example \
 		-p ${TEPIC_DIR}/PWMs/1.0/pwm_vertebrates_jaspar_uniprobe_original.PSEM \
 		-a ${TEPIC_DIR}/Test/example_annotation.gtf \
 		-w 3000 \
 		-e FALSE

# TEPIC INVOKE pipeline
# =====config==========
define TEPIC_INVOKE_CFG
#The following parameters should be changed by the user
open_regions="${INVOKE_DIR}/ExampleData/S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1.bed"
#Note that the brackets at the beginning and at the end are compulsory.
gene_Expression_Data=("${INVOKE_DIR}/ExampleData/S001S712_gene_quantification_rsem_grape2_crg_GRCh38_20150622_chr1.txt")
outputDirectory="${INVOKE_OUT}"
preComputedTEPIC="${INVOKE_DIR}/TestRun/Affinities/S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1_TEPIC_06_08_17_15_34_52_Decay_Affinity_Gene_View_Filtered.txt"

#Parameters of TEPIC
path=${TEPIC_DIR}
referenceGenome="${INVOKE_DIR}/ExampleData/Homo_sapiens.GRCh38.dna.chromosome.1.fa"
pwm=$${path}/PWMs/1.0/pwm_vertebrates_jaspar_uniprobe_original.PSEM
cores_TEPIC="${NCPU}"
geneAnnotation=""
window="10000"
peakFeature="FALSE"
decay="TRUE"
chrPrefix="FALSE"
#Note that the brackets at the beginning and at the end are compulsory.
coverageFile=("")
coverageColumn=("")

#Parameters of the linear regression
scriptPath="${INVOKE_DIR}/Scripts"
cores_R="${NCPU}"
alpha_Step_Size="0.1"
regularisation="E"
innerCV="6"
outerCV="3"
testsize="0.2"
performance="TRUE"
randomise="FALSE"
fTest="FALSE"
endef
export TEPIC_INVOKE_CFG := $(TEPIC_INVOKE_CFG)
# =====config==========

.PHONY: tepic_invoke
tepic_invoke: ## run TEPIC INVOKE pipeline on it's example data
	@echo "Running TEPIC INVOKE pipeline"
	mkdir -p ${INVOKE_OUT}
	echo "$${TEPIC_INVOKE_CFG}" > ${INVOKE_OUT}/invokeAnalysis.cfg
	${INVOKE_BIN} ${INVOKE_OUT}/invokeAnalysis.cfg

# TEPIC DYNAMITE pipeline
# =====config==========
define TEPIC_DYNAMITE_CFG
#The following parameters should be changed by the user
open_regions_Group1="${DYNAMITE_DIR}/ExampleData/Hf03_BlEM_NOMe.bed"
open_regions_Group2="${DYNAMITE_DIR}/ExampleData/Hf03_BlTN_NOMe.bed"
differential_Gene_Expression_Data="${DYNAMITE_DIR}/ExampleData/TEM_vs_TN0_000001.txt"
outputDirectory="${DYNAMITE_OUT}"

#To reuse existing TEPIC results on two different groups specfiy the path to the group1 and group2 folders here. Note that mean affinities and ratios are recomputed
existing_TEPIC_Results_Group1="${DYNAMITE_DIR}/TestRun/Affinities/group1/"
existing_TEPIC_Results_Group2="${DYNAMITE_DIR}/TestRun/Affinities/group2/"

#Parameters of TEPIC
path=${TEPIC_DIR}
referenceGenome=""
pwm="${DYNAMITE_DIR}/ExampleData/T-Cell-Example.PSEM"
cores_TEPIC=${NCPU}
geneAnnotation=""
window="50000"
coverage_Files_Group1=("")
coverage_Files_Group2=("")
coverage_Columns_Group1=("")
coverage_Columns_Group2=("")
decay="TRUE"
peakFeatures="TRUE"
chrPrefix="FALSE"

#Parameters of the logistic regression
scriptPath="${DYNAMITE_DIR}/Scripts"
outerCV=3
innerCV=6
alpha_Step_Size=0.1
coresR=${NCPU}
performance="TRUE"
randomise="FALSE"
endef
export TEPIC_DYNAMITE_CFG := $(TEPIC_DYNAMITE_CFG)
# =====config==========

.PHONY: tepic_dynamite
tepic_dynamite: ## run TEPIC DYNAMITE pipeline on it's example data
	@echo "Running TEPIC DYNAMITE pipeline"
	mkdir -p ${DYNAMITE_OUT}
	echo "$${TEPIC_DYNAMITE_CFG}" > ${DYNAMITE_OUT}/DYNAMITE.cfg
	${DYNAMITE_BIN} ${DYNAMITE_OUT}/DYNAMITE.cfg
	Rscript ${DYNAMITE_DIR}/Scripts/generateFeaturePlots.R ${DYNAMITE_OUT}/ Peak_Counts TEM TN

# Based on http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.DEFAULT_GOAL := help
.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
