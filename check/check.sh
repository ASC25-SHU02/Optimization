#!/bin/bash
echo ""

function clear_files() {
	rm correlation.tsv  detected.tsv
}
clear_files

python utils/get_intersection.py

bash ./utils/precision.sh

python utils/get_correlation.py

bash ./utils/correlation.sh
