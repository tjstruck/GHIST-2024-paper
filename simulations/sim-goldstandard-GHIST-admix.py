#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-admix
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


challenge_name = "GHIST-admixture"
os.makedirs(challenge_name, exist_ok=True)

def msprime_admixture(Na, p):
    (nuExt, nuAnc, nuMod1_1, nuMod1_2, nuMod2_1, nuMod2_2, TSplitAnc, TSplitMod, TAdmix1, TAdmix2, m) = p
    dem = msprime.Demography()
    dem.add_population(name="common_ancestor", initial_size=Na, initially_active=True)
    dem.add_population(name="Modern_pop1", initial_size=nuMod1_2)
    dem.add_population(name="Modern_pop2", initial_size=nuMod2_2)
    dem.add_population(name="AncOut1", initial_size=nuExt, default_sampling_time=4000)
    dem.add_population(name="AncOut2", initial_size=nuExt, default_sampling_time=5000)
    dem.add_population(name="AncPop1", initial_size=nuExt, default_sampling_time=1000)#5200)
    dem.add_population(name="AncOut3", initial_size=nuExt, default_sampling_time=1000)
    dem.add_population(name="AncPop2", initial_size=nuExt, default_sampling_time=700)
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population(name="ancient_population", initial_size=nuAnc)

    dem.add_population_parameters_change(time=TSplitAnc*1.5, population="ancestral", initial_size=Na*.5)
    dem.add_population_parameters_change(time=TSplitAnc*1.2, population="ancestral", initial_size=Na*1.5)
    dem.add_population_parameters_change(time=TSplitAnc, population="ancestral", initial_size=Na*2)

    dem.add_population_split(time=TSplitAnc, derived=["common_ancestor", "ancient_population"], ancestral="ancestral")

    dem.add_population_split(time=TSplitMod, derived=["Modern_pop1", "Modern_pop2"], ancestral="ancient_population")

    dem.add_population_split(time=13000, derived=["AncOut1"], ancestral="common_ancestor")
    dem.add_population_split(time=11000, derived=["AncOut2"], ancestral="common_ancestor")
    dem.add_population_split(time=8000, derived=["AncPop1"], ancestral="common_ancestor")
    dem.add_population_split(time=10000, derived=["AncOut3"], ancestral="common_ancestor")
    dem.add_population_split(time=7000, derived=["AncPop2"], ancestral="common_ancestor")


    dem.add_mass_migration(time=TAdmix1, source="Modern_pop1", dest="AncPop1", proportion=pulse_size1)
    dem.add_mass_migration(time=TAdmix2, source="Modern_pop2", dest="AncPop2", proportion=pulse_size2)

    dem.add_migration_rate_change(time=0, rate=m, source="Modern_pop1", dest="Modern_pop2")
    dem.add_migration_rate_change(time=0, rate=m, source="Modern_pop2", dest="Modern_pop1")
    dem.add_migration_rate_change(time=0, rate=m, source="AncPop1", dest="AncPop2")
    dem.add_migration_rate_change(time=0, rate=m, source="AncPop2", dest="AncPop1")

    dem.add_migration_rate_change(time=0, rate=m, source="AncOut3", dest="AncPop2")
    dem.add_migration_rate_change(time=0, rate=m, source="AncPop2", dest="AncOut3")

    dem.add_migration_rate_change(time=0, rate=m, source="AncOut1", dest="AncPop1")
    dem.add_migration_rate_change(time=0, rate=m, source="AncPop1", dest="AncOut1")

    dem.add_migration_rate_change(time=0, rate=m, source="AncPop1", dest="AncOut2")
    dem.add_migration_rate_change(time=0, rate=m, source="AncOut2", dest="AncPop1")

    dem.add_migration_rate_change(time=0, rate=m, source="AncOut1", dest="AncOut2")
    dem.add_migration_rate_change(time=0, rate=m, source="AncOut2", dest="AncOut1")
    dem.add_migration_rate_change(time=0, rate=m, source="AncOut1", dest="AncOut3")
    dem.add_migration_rate_change(time=0, rate=m, source="AncOut3", dest="AncOut1")
    dem.add_migration_rate_change(time=0, rate=m, source="AncOut2", dest="AncOut3")
    dem.add_migration_rate_change(time=0, rate=m, source="AncOut3", dest="AncOut2")

    dem.add_population_parameters_change(time=TAdmix1*.8, population="Modern_pop1", initial_size=nuMod1_1)
    dem.add_population_parameters_change(time=TAdmix2*.8, population="Modern_pop2", initial_size=nuMod2_1)
    return dem

Na = 32671 # ancestral pop size
nuExt = 13249 # extinc pop size
nuAnc = 5083 # common ancestor pop size

nuMod1_1 = 2231.0 # modern pop 1 post-split size
nuMod1_2 = 9025.0 # modern pop 1 current size

nuMod2_1 = 1293.0 # modern pop 2 post-split size
nuMod2_2 = 6962.0 # modern pop 2 current size

TSplitAnc = 15090.0 # generations ago for split of common ancestor and extinct population
TSplitMod = 1758.0 # generations ago for split of modern populations

TAdmix1 = 1566 # generations ago for admixture into modern pop 1
TAdmix2 = 883 # generations ago for admixture into modern pop 2

pulse_size1 = 0.011
pulse_size2 = 0.002

m = 3.14e-05

p = (nuExt, nuAnc, nuMod1_1, nuMod1_2, nuMod2_1, nuMod2_2, TSplitAnc, TSplitMod, TAdmix1, TAdmix2, m)

ns = {"Modern_pop1":20, "Modern_pop2":16, "AncOut1":3, "AncOut2":2, "AncPop1": 2, "AncOut3":2, "AncPop2": 1} # individuals
ploidy = 2 # diploid
mut = 1.29e-8 # mutation rate

recomb = 1.38e-08
# msprime demography model with dadi parameters
dem = msprime_admixture(Na, p)
dem.sort_events()

ts = msprime.sim_ancestry(
    samples=ns, 
    demography=dem, 
    sequence_length=2.5e8, 
    recombination_rate=recomb, 
    ploidy=ploidy,
    random_seed=12
    )

mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=False, random_seed=5566)

vcf_fi = open(f"{challenge_name}/{challenge_name}.vcf","w")
vcf_fi.write(mts.as_vcf())
vcf_fi.close()

fi = open(f"{challenge_name}/{challenge_name}.popfile","w")

for i in range(1,ns['Modern_pop1']+1):
    fi.write(f"Modern1.{i}\tModern1\n")

for i in range(1,ns['Modern_pop2']+1):
    fi.write(f"Modern2.{i}\tModern2\n")

for i in range(1,ns['AncOut1']+1):
    fi.write(f"AncSite1.{i}\tAnc1\n")

for i in range(1,ns['AncOut2']+1):
    fi.write(f"AncSite2.{i}\tAnc2\n")

for i in range(1,ns['AncPop1']+1):
    fi.write(f"AncSite3.{i}\tAnc3\n")

for i in range(1,ns['AncOut3']+1):
    fi.write(f"AncSite4.{i}\tAnc4\n")

for i in range(1,ns['AncPop2']+1):
    fi.write(f"AncSite5.{i}\tAnc5\n")

fi.close()