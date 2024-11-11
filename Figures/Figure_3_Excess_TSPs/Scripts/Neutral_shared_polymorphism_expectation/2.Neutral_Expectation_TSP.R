# Neutral expectation of number of TSPs (upper bound)
dAB=0.0535 # Between species dxy
dA=0.0186 # within NAm. D. pulex dxy
dB=0.00413 # within Euro D. pulex dxy
num.snps=225734 # Total number of BUSCO Gene SNPs

e1=exp(-((dAB-max(dA,dB))/dA))
e2=exp(-((dAB-max(dA,dB))/dB))
prop.tsp=e1*e2

# Expected Number of TSPs
(prop.tsp*num.snps)

# Actual Number of TSPs - Synonymous
emp=4075

# Enrichment
log2(emp/(prop.tsp*num.snps))

# Prob. of TSP coalescent
split=10000000
ne=700000
exp(-split/(2*ne))

0.0351*637739
