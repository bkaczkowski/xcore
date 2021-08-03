PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
IMAGE_BIN = $(PROJECT_DIR)/IMAGE/IMAGE.pl
NCPU = 6

.PHONY: install
install: ## download programs and set up environment
	@echo "setting up IMAGE"
	wget -c http://bioinformatik.sdu.dk/solexa/webshare/IMAGE/IMAGE_v1.1.tar.gz
	tar -xzf IMAGE_v1.1.tar.gz
	wget -c -P IMAGE/ https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
	gunzip IMAGE/mm9.fa.gz

.PHONY: IMAGE
image: ## run IMAGE on it's example data
	@echo "Running IMAGE"
	${IMAGE_BIN} \
		-region IMAGE/examples/Enhancers.txt \
		-expression IMAGE/examples/GeneExpression.txt \
		-fasta IMAGE/mm9.fa \
		-RNADesign 1 1 2 2 3 3 4 4 \
		-EnhancerDesign 1 1 2 2 3 3 4 4 \
		-p ${NCPU} \
		-n IMAGE_test
	mv IMAGE_test.R IMAGE_test.Rdata

# Based on http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.DEFAULT_GOAL := help
.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
