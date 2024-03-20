grep "CV" *[0-9].out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > daphnia.cv.error
grep "CV" *noobtusa.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > daphnia.noobtusa.cv.error
grep "CV" *sub.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > daphnia.sub.cv.error

# Go through each population
for Individuals_to_keep in /scratch/csm6hg/daphnia_phylo/admixture/popfile.Daphnia*; do

  # Progress message
  echo $Individuals_to_keep

  # Output name
  spp=$( echo $Individuals_to_keep | cut -f6 -d"/" | cut -f2-4 -d"." )
  echo $spp

  # For each independent species
  grep "CV" *${spp}.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > ${spp}.sub.cv.error

# Finish
done
