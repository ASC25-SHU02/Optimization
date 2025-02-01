#!/bin/bash
echo ""
python utils/get_intersection.py

bash ./utils/precision.sh

python utils/get_intersection.py

bash ./utils/correlation.sh