#!/bin/bash

python utils/get_intersection.py

sh ./utils/precision.sh

python utils/get_intersection.py

sh ./utils/correlation.sh