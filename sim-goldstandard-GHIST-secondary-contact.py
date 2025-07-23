#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-sec-cont
#SBATCH --output=hpc_output/%x-%A.out
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjstruck@arizona.edu
###SBATCH --partition=windfall
#SBATCH --partition=high_priority
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --cpus-per-task=1
###SBATCH --mem-per-cpu=5gb
###SBATCH --constraint=hi_mem
###SBATCH --mem-per-cpu=32gb
#SBATCH --time=72:00:00
###SBATCH --array=1-3

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

challenge_name = "GHIST-secondary-contact"
os.makedirs(challenge_name, exist_ok=True)

def msprime_split_isolation_then_migtation(Na, p):
    (nu1, nu2, nu2Cont, TAncSplit, TCont, m) = p
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population(name="mainland", initial_size=Na*nu1)
    dem.add_population(name="island", initial_size=Na*nu2)

    dem.add_population_parameters_change(time=TAncSplit*1.2, population="ancestral", initial_size=Na*.5)
    dem.add_population_parameters_change(time=TAncSplit*1.1, population="ancestral", initial_size=Na*1.5)
    dem.add_population_parameters_change(time=TAncSplit*1.05, population="ancestral", initial_size=Na*2)

    dem.add_population_split(time=TAncSplit, derived=["mainland", "island"], ancestral="ancestral")
    dem.add_population_parameters_change(time=TAncSplit*0.8, population="mainland", initial_size=Na*0.8)
    dem.add_population_parameters_change(time=TAncSplit*0.5, population="mainland", initial_size=Na*0.5)
    dem.add_population_parameters_change(time=0, population="mainland", initial_size=Na*nu1)

    # dem.set_symmetric_migration_rate(["mainland", "island"], 0)
    dem.add_migration_rate_change(time=TCont, rate=0, source="mainland", dest="island")
    dem.add_migration_rate_change(time=TCont, rate=0, source="island", dest="mainland")
    dem.add_migration_rate_change(time=0, rate=m, source="mainland", dest="island")
    dem.add_migration_rate_change(time=0, rate=m, source="island", dest="mainland")

    dem.add_population_parameters_change(time=TCont, population="island", initial_size=Na*nu2)
    dem.add_population_parameters_change(time=TCont*0.45, population="island", initial_size=Na*nu2*1.5)
    dem.add_population_parameters_change(time=0, population="island", initial_size=Na*nu2Cont)
    return dem


Na = 120000 # ancestral pop size
nu1 = 2
nu2 = 0.1
nu2Cont = 0.3
TAncSplit = 23000
TCont = 1277 
m = 5e-5
p = (nu1, nu2, nu2Cont, TAncSplit, TCont, m)

ns = {"mainland":22, "island":8} # individuals
ploidy = 2 # diploid
# based on https://github.com/popsim-consortium/stdpopsim/blob/main/stdpopsim/catalog/MusMus/species.py
mut = 5.4e-9 # mutation rate
# mean recomb from stdpopsim
recomb = 5.386e-09
# msprime demography model with dadi parameters
dem = msprime_split_isolation_then_migtation(Na, p)
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

vcf_fi = open(f"{challenge_name}/{challenge_name}{less}.vcf","w")
vcf_fi.write(mts.as_vcf())
vcf_fi.close()
