import moments
import numpy as np
import yaml
import pickle
import sys
import os

# 2DSFS inference

# Import config YAML file into global dictionary.
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

def model_func(params, ns, pop_ids=None):
    if modeli == "split_mig":
        return moments.Demographics2D.split_mig(params, ns, pop_ids=pop_ids).fold()
    elif modeli == "split_no_mig":
        nu1, nu2, T = params
        sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
        fs = moments.Spectrum(sts)
        fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
        fs.integrate([nu1, nu2], T, m=np.zeros([2, 2]))
        fs.pop_ids = pop_ids
        return fs.fold()
    else:
        raise "Model name not recognized."

# Load VCF as serialized data_dict object.
with open(yd["data_dict_pickle_file"], "rb") as f:
    data_dict = pickle.load(f)
# Convert data_dict object into Spectrum object.

# Metadata from command line
proj = sys.argv[1]
outNamei = sys.argv[2]
modeli = sys.argv[3]

# Run moments
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=yd["pop_names"],
                                     projections=[int(proj)] * 2,
                                     polarized=False)

print("SFS loaded.")

# Setup for moments inference.
model_to_num_of_params = {"split_mig": 4, "split_no_mig": 3, "split_mig_asym": 5}
lower_bound = [1e-4 for i in range(model_to_num_of_params[modeli])]
upper_bound = [10 for i in range(model_to_num_of_params[modeli])]

ns = [i - 1 for i in fs.shape]

# Perform inference "rep" many times and choose the parameter estimates from the
# inference run with the greatest likelihood.
for rep in range(yd["num_of_inference_repeats"]):
    print("Rep:", rep)
    # Randomly generate initial guess for parameter estimates
    params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(len(lower_bound))]
    # Perform inference to get parameter estimates and assess model performance
    popt = moments.Inference.optimize_log(params, fs, model_func,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound,
                                          verbose=yd["verbosity"])
    model = model_func(popt, ns)
    theta = moments.Inference.optimal_sfs_scaling(model, fs)
    log_likelihood = moments.Inference.ll_multinom(model, fs)

    # Convert from coalescent units.
    conversion_coeff = theta / (4 * yd["mutation_rate"] * yd["seq_len"])
    popt[0] *= conversion_coeff
    popt[1] *= conversion_coeff
    popt[2] *= 2 * conversion_coeff
    if modeli == "split_mig":
        popt[3] /= 2 * conversion_coeff

    # Print to output file.
    with open(yd["output_file"] + "." + outNamei, "a") as f:
        f.write(modeli + "\t")
        f.write(proj + "\t")
        f.write(str(theta) + "\t")
        f.write(str(log_likelihood) + "\t")
        for param in popt:
            f.write(str(param) + "\t")
        if modeli == "split_no_mig":
            f.write("NA")
        f.write("\n")


