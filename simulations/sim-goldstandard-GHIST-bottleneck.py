#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-bottleneck
#SBATCH --output=hpc_output/%x-%A.out
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjstruck@arizona.edu
###SBATCH --partition=windfall
#SBATCH --partition=high_priority
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
###SBATCH --mem-per-cpu=5gb
###SBATCH --constraint=hi_mem
###SBATCH --mem-per-cpu=32gb
#SBATCH --time=1:00:00
###SBATCH --array=1-50

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


challenge_name = "GHIST-bottleneck"
os.makedirs(challenge_name, exist_ok=True)

seq_l = "1e8"

# if seq_l != "1e8":
seq_tag = "."+seq_l
# else:
seq_tag = ''

def msprime_wisent_bottleneck(Na, p):
    (nu, T) = p
    dem = msprime.Demography()
    dem.add_population(name="wisent", initial_size=Na*nu)
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population_split(time=2*Na*T, derived=["wisent"], ancestral="ancestral")
    return dem

nu = 0.08
T = 0.01

p = (nu, T)

Na = 14182
ns = {"wisent":20} # individuals
ploidy = 2 # diploid
mut = 1.26e-8 # mutation rate

# 1e-8 recombination per base seems roughtly right 0.0099±0.0052 and 0.0088±0.0053 per 1MB in two species of cattle
recomb = 1.007e-08
# msprime demography model with dadi parameters
dem = msprime_wisent_bottleneck(Na, p)
dem.sort_events()

os.makedirs(f"{challenge_name}", exist_ok=True)

ts = msprime.sim_ancestry(
    samples=ns, 
    demography=dem, 
    sequence_length=1e8, 
    recombination_rate=recomb, 
    ploidy=ploidy,
    random_seed=12
    )

mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=False, random_seed=5566)
mts.num_mutations

vcf_fi = open(f"{challenge_name}{seq_tag}/{challenge_name}.vcf","w")
vcf_fi.write(mts.as_vcf())
vcf_fi.close()

fi = open(f"{challenge_name}/{challenge_name}.popfile","w")
for i in range(1,21):
    fi.write(f"wisent_{i}\twisent\n")
fi.close()