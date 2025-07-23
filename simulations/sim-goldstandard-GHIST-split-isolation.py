#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-split-iso
#SBATCH --output=hpc_output/%x-%A.out
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjstruck@arizona.edu
###SBATCH --partition=windfall
#SBATCH --partition=high_priority
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
###SBATCH --mem-per-cpu=5gb
###SBATCH --constraint=hi_mem
###SBATCH --mem-per-cpu=32gb
#SBATCH --time=72:00:00
###SBATCH --array=1-4

import dadi
import msprime
import demes
import random
import ast
import sys,os
# Python-specific stuff
print('Script running\n')
if 'SLURM_SUBMIT_DIR' in os.environ:
#       # Set module search path so it will search in qsub directory
      sys.path.insert(0, os.environ['SLURM_SUBMIT_DIR'])
#       # Set current working directory to qsub directory
      # os.chdir(os.environ['SLURM_SUBMIT_DIR'])
# Which process am I?
process_ii = int(os.environ.get('SLURM_ARRAY_TASK_ID',1))-1
print(process_ii)


challenge_name = "GHIST-split-isolation"
os.makedirs(challenge_name, exist_ok=True)

def msprime_split_isolation(Na, p):
    (nu1, nu2, T) = p
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population(name="east", initial_size=Na*nu1)
    dem.add_population(name="west", initial_size=Na*nu2)
    dem.add_population_split(time=T, derived=["east", "west"], ancestral="ancestral")

    return dem

Na = 100000 # ancestral pop size

nu1 = 1.3
nu2 = 0.2#random.random()*10
T = 13333
p = (nu1, nu2, T)


ns = {"east":22, "west":18} # individuals
ploidy = 2 # diploid
# based on https://github.com/popsim-consortium/stdpopsim/blob/main/stdpopsim/catalog/MusMus/species.py
mut = 5.7e-9 # mutation rate
# mean recomb from stdpopsim
recomb = 5.386e-09
# msprime demography model with dadi parameters
dem = msprime_split_isolation(Na, p)
dem.sort_events()


ts = msprime.sim_ancestry(
    samples=ns, 
    demography=dem, 
    sequence_length=1e8, 
    recombination_rate=recomb, 
    ploidy=ploidy,
    random_seed=12
    )
ts

mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=False, random_seed=5566)
mts.num_mutations

vcf_fi = open(f"{challenge_name}/{challenge_name}.vcf","w")
vcf_fi.write(mts.as_vcf())
vcf_fi.close()

fi = open(f"{challenge_name}/{challenge_name}.popfile","w")

for i in range(1,ns['east']+1):
    fi.write(f"east_{i}\teast\n")
for ii in range(1, ns['west']+1):
    fi.write(f"west_{ii}\twest\n")
fi.close()