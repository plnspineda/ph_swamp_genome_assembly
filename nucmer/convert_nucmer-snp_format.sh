#!/bin/bash
### PSPINEDA 2023
### for converting the SNP alignment from nucmer (from this run: `show-snps -Clr "$OUTFILE" > "$FINOUT"`) to a machine readable file

input=$1
base=$(basename $input)
tail -n+6 "$input" | sed 's/|//g' | awk '{$1=$1}1' | sed 's/ /\t/g' > tmp && mv tmp "$base".bed
