#!/bin/bash


#
#SBATCH --job-name=outSchProxy
#SBATCH --output=./SlurmOutput/outSchProxy-%j.out
#
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --exclude=node[120-128]



#julia outageSchedulingProxy.jl ./data/IEEE-RTS96 ./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData.csv 16 1 1 1 10 1 182 1


julia outageSchedulingProxy.jl ./data/IEEE-RTS96 ./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData.csv 16 1 1 1 32 1 182 2
