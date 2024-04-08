#!/bin/bash

# Input BED file
input_tab_delimited_file="nampulex_chroms.bed"

# Output BED file
output_bed_file="ChrEnd500bps.nampulex.bed"

# Number of base pairs to extract from the beginning and end of each chromosome
bp_to_extract=500

# Create a temporary directory to store intermediate files
temp_dir=$(mktemp -d)

# Use awk to extract the first and last 500 bp of each chromosome in tab-delimited format
awk -v bp_to_extract="$bp_to_extract" 'BEGIN {OFS="\t"} {
    if (!($1 in chrom_start)) {
        chrom_start[$1] = $2;
    }
    chrom_end[$1] = $3;
}
END {
    for (chrom in chrom_start) {
        # Extract the first 500 bp
        print chrom, chrom_start[chrom], chrom_start[chrom] + bp_to_extract;
        # Extract the last 500 bp
        print chrom, chrom_end[chrom] - bp_to_extract, chrom_end[chrom];
    }
}' "$input_bed_file" > "$temp_dir/extracted.tab"

# Sort the extracted TAB-delimited file by chromosome and position
sort -k1,1 -k2,2n "$temp_dir/extracted.tab" > "$output_bed_file"

# Clean up the temporary directory
rm -r "$temp_dir"

echo "Extracted regions saved to $output_tab_delimited_file"

