#!/bin/bash
echo "=====Correlation====="
awk '{
sumXY += $1 * $2;
sumX += $1;
sumY += $2;
sumX2 += $1 * $1;
sumY2 += $2 * $2;
n++;
} END {
print (n * sumXY - sumX * sumY) / (sqrt(n * sumX2 - sumX^2) * sqrt(n * sumY2 - sumY^2))
}' correlation.tsv

echo "====================="
