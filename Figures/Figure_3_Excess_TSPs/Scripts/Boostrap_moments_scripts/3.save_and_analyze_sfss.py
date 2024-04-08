# Packages to import
import moments
import matplotlib
import numpy as np
import pickle
from math import log
import sys
import os

# Given moments' parameter estimates, reconstruct SFSs corresponding to optimal models.
# Serialize the SFSs as .npy files and save them as CSVs for later use in R visualization.
# Calculate summary statistics like shared allele proportion (SAP) and normalized log-likelihood (ll)
# and save them to an output file.

# Setup and parameters
dir = "/scratch/csm6hg/moments/"
vcf = "filt"
pickle_file = dir + "daphnia.filt.mlg.genome.11.18.22_data_dict.pickle"
csv_dir = dir + "sfss/sfs_csvs/"
sfs_dir = dir + "sfss/sfs_npys/"
pop_ids = ["Daphnia.pulex.NorthAmerica", "Daphnia.pulex.Europe"]
prior_estimated_params = [7e5, 2e5, 1e7 / 2, 1e-8 * 2]

# Old optimized values
# split_mig_params = [6.758185, 1.1126711, 9.9724427, 0.008836859]
# split_no_mig_params = [1.127693, 0.2289919, 0.7780228, 0]

# Model
modeli = sys.argv[1]
print(modeli)

# Optimized parameters
OptParam = sys.argv[2]
# OptParam="6.758185,1.1126711,9.9724427,0.008836859"
OptParam = tuple(OptParam.split(','))
OptParam = tuple([float(num) for num in OptParam])
print(OptParam)

# Replicate
repli = sys.argv[3]
print(repli)

# Output file
output_file = dir + "output/sfs_statistics" + "." + modeli + "." + repli + ".txt"

# modeli="split_mig": OptParam=6.758185,1.1126711,9.9724427,0.008836859: repli=3

# Summary Statistic Functions
def save_sfs_as_csv_and_npy(sfs, sfs_name):
    # Save SFS as CSV
    np.savetxt(csv_dir + sfs_name + ".csv", sfs, delimiter=",")
    # Save SFS as numpy binary
    with open(sfs_dir + sfs_name + ".npy", "wb") as fout:
        np.save(fout, sfs.data)

def get_AIC(sfs, sfs_empirical, k):
    # k is the number of model parameters.
    ll = moments.Inference.ll_multinom(sfs, sfs_empirical)
    AIC = k * 2 - 2 * ll
    return AIC

def get_BIC(sfs, sfs_empirical, k, sfs_size):
    # k is the number of model parameters.
    # The number of elements in a square folded SFS with the absent allele corner
    # masked is the ith triangular number minus 1, where i is the sample size.
    n = sfs_size * (sfs_size + 1) / 2 - 1
    ll = moments.Inference.ll_multinom(sfs, sfs_empirical)
    BIC = k * log(n) - 2 * ll
    return BIC

def get_normalized_ll(sfs, sfs_empirical):
    ll = moments.Inference.ll_multinom(sfs, sfs_empirical / np.sum(sfs_empirical))
    return ll

def get_shared_allele_prop(sfs, maf_threshold=0.01):
    # Iterate through sfs, summing element values, but skipping entries corresponding
    # to alleles that are not sufficiently shared according to "maf_threshold"
    shared_allele_count = 0
    coordinate_thresholds = np.array([(len - 1) * maf_threshold + 1 for len in sfs.shape]).astype(int)

    # Create NumPy iterator object
    it = np.nditer(sfs, flags=['multi_index'])
    for i in it:
        shared_allele = True
        # Check for whether element should be skipped because it fails to cross
        # MAF threshold
        for j, coord in enumerate(it.multi_index):
            if coord < coordinate_thresholds[j]:
                shared_allele = False
                break
        # If it crosses the MAF threshold to be considered "shared", add it to the counter
            if shared_allele:
                shared_allele_count += i
                shared_allele_prop = shared_allele_count / np.sum(sfs)
                return shared_allele_prop

def write_output(fout, outputs):
    for output in outputs:
        fout.write(str(output) + "\t")
    fout.write("\n")

# Load serialized data_dict object, which contains the same information as the VCF
# file, but formatted for moments
with open(pickle_file, "rb") as f:
    data_dict = pickle.load(f)

# ns gives dimensions of SFS as [ns[0] + 1, ns[1] + 1]
with open(output_file, "a") as fout:
    # Write header
    #fout.write("n\tsfs\tAIC\tBIC\tSAP\n")

    # Go through various projection sizes
    for ns in [[20, 20], [100, 100]]:

        # Empirical
        sfs_empirical = moments.Spectrum.from_data_dict(data_dict, pop_ids=pop_ids,
                                                  projections=ns,
                                                  polarized=False)
        save_sfs_as_csv_and_npy(sfs_empirical, "sfs_empirical_" + str(ns[0]) + "." + repli)
        outputs = [vcf, ns[0], "sfs_empirical", "NA", "NA", "NA", get_shared_allele_prop(sfs_empirical)]
        write_output(fout, outputs)

        # Run Prior
        sfs_from_ests = moments.Demographics2D.split_mig(prior_estimated_params, ns,
                                                   pop_ids=pop_ids).fold()
        sfs_from_ests *= moments.Inference.optimal_sfs_scaling(sfs_from_ests, sfs_empirical)
        save_sfs_as_csv_and_npy(sfs_from_ests, "sfs_from_ests_" + str(ns[0]) + "." + repli)
        outputs = [vcf, ns[0], "sfs_from_ests",
                   get_AIC(sfs_from_ests, sfs_empirical, 4),
                   get_BIC(sfs_from_ests, sfs_empirical, 4, ns[0] + 1),
                   get_normalized_ll(sfs_from_ests, sfs_empirical),
                   get_shared_allele_prop(sfs_from_ests)]
        write_output(fout, outputs)

        # Conditional run "Split + Mig" model
        if modeli == "split_mig":
            print("Split_mig model running")
            sfs_split_mig_model = moments.Demographics2D.split_mig(OptParam, ns, pop_ids=pop_ids).fold()
            sfs_split_mig_model *= moments.Inference.optimal_sfs_scaling(sfs_split_mig_model, sfs_empirical)
            save_sfs_as_csv_and_npy(sfs_split_mig_model, f"sfs_split_mig_model_{ns[0]}.{repli}")
            outputs = [vcf, ns[0], "sfs_split_mig_model",
                    get_AIC(sfs_split_mig_model, sfs_empirical, 4),
                    get_BIC(sfs_split_mig_model, sfs_empirical, 4, ns[0] + 1),
                    get_normalized_ll(sfs_split_mig_model, sfs_empirical),
                    get_shared_allele_prop(sfs_split_mig_model)]
            write_output(fout, outputs)

        # Conditional run "Split + No Mig" model
        elif modeli == "split_no_mig":
            print("Split_no_mig model running")
            sfs_split_no_mig_model = moments.Demographics2D.split_mig(OptParam, ns, pop_ids=pop_ids).fold()
            sfs_split_no_mig_model *= moments.Inference.optimal_sfs_scaling(sfs_split_no_mig_model, sfs_empirical)
            save_sfs_as_csv_and_npy(sfs_split_no_mig_model, f"sfs_split_no_mig_model_{ns[0]}.{repli}")
            outputs = [vcf, ns[0], "sfs_split_no_mig_model",
                    get_AIC(sfs_split_no_mig_model, sfs_empirical, 3),
                    get_BIC(sfs_split_no_mig_model, sfs_empirical, 3, ns[0] + 1),
                    get_normalized_ll(sfs_split_no_mig_model, sfs_empirical),
                    get_shared_allele_prop(sfs_split_no_mig_model)]
            write_output(fout, outputs)