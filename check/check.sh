#!/bin/bash
echo ""

function clear_files() {
	if [ -f "$1" ]; then  
  		rm "$1" 
	fi 
}
clear_files correlation.tsv
clear_files detected.tsv

python utils/get_intersection.py

bash ./utils/precision.sh

python utils/get_correlation.py

bash ./utils/correlation.sh
