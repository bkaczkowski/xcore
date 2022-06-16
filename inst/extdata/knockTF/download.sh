#!/bin/bash
cat encode_datasets.txt | xargs -n 1 -P 6 wget -c
