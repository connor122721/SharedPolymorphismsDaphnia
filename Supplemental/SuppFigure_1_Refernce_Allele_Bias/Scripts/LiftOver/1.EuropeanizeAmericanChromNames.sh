#!/usr/bin/env bash

# We are changing the north american chromosome names to their corresponding chromosome in euro pulex
cd /scratch/csm6hg/data

# Rename chromosomes in NAM pulex
sed -i 's/NC_060026.1/Scaffold_2217_HRSCAF_2652/g' american_to_european_chredit.liftOver
sed -i 's/NC_060023.1/Scaffold_7757_HRSCAF_8726/g' american_to_european_chredit.liftOver
sed -i 's/NC_060019.1/Scaffold_9199_HRSCAF_10755/g' american_to_european_chredit.liftOver
sed -i 's/NC_060018.1/Scaffold_9198_HRSCAF_10754/g' american_to_european_chredit.liftOver
sed -i 's/NC_060024.1/Scaffold_6786_HRSCAF_7541/g' american_to_european_chredit.liftOver
sed -i 's/NC_060021.1/Scaffold_9200_HRSCAF_10757/g' american_to_european_chredit.liftOver
sed -i 's/NC_060028.1/Scaffold_2158_HRSCAF_2565/g' american_to_european_chredit.liftOver
sed -i 's/NC_060017.1/Scaffold_1931_HRSCAF_2197/g' american_to_european_chredit.liftOver
sed -i 's/NC_060020.1/Scaffold_9197_HRSCAF_10753/g' american_to_european_chredit.liftOver
sed -i 's/NC_060025.1/Scaffold_1863_HRSCAF_2081/g' american_to_european_chredit.liftOver
sed -i 's/NC_060027.1/Scaffold_9201_HRSCAF_10758/g' american_to_european_chredit.liftOver
sed -i 's/NC_060022.1/Scaffold_2373_HRSCAF_2879/g' american_to_european_chredit.liftOver
