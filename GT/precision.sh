# Step 1: Find rows in detected.tsv that match rows in true.tsv based on the first three columns
matches=$(awk 'NR==FNR {a[$1,$2,$3]=1; next} ($1,$2,$3) in a' "true.tsv" "detected.tsv")

# Step 2: Count the number of matching rows (true positives)
true_positives=$(echo "$matches" | wc -l)

# Step 3: Count the total number of rows in detected.tsv
total_detected=$(wc -l < "detected.tsv")

# Step 4: Calculate precision as a percentage
precision=$(awk -v total="$total_detected" -v tp="$true_positives" 'BEGIN {printf "%.2f", (tp/total)*100}')

# Output the precision
echo "Precision: $precision%"
