# Modules
# module load anaconda/2023.07-py3.11; conda activate msprime_env; python

# Libraries
import msprime
import seaborn as sns
import matplotlib.pyplot as plt

# Model ancestry
ts = msprime.sim_ancestry(5,
                          population_size=10000,
                          sequence_length=130000000,
                          recombination_rate=5e-9,
                          random_seed=1235123515)

# Making lists of time-to the most recent common ancestors
tmrca = []
span = []
for t in ts.trees():
    for r in t.roots:
        tmrca.append(ts.node(r).time / 10000)
        span.append(t.span)

# Plotting TMRCA distributions
fig = sns.histplot(tmrca)
fig.set_xlabel("Scaled TMRCA")
plt.savefig("tmrca.pdf")
plt.clf()

fig = sns.scatterplot(x=tmrca, y=span)
fig.set_xlabel("Scaled TMRCA")
fig.set_ylabel("Tree span (bp)")
plt.savefig("tmrca_vs_span.pdf")