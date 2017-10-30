#!/bin/bash


#
#SBATCH --job-name=genMicro
#SBATCH --output=./SlurmOutput/genMicro-%j.out
#
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --exclude=node[120-128]

#julia script, called with: spec_filename ntasks cpus-per-task nodes affected_branch first_micro last_micro first_day last_day
julia generateMicroscenarios.jl ./data/IEEE-RTS96 16 1 1 65 128 1 182
