#!/bin/bash


#
#SBATCH --job-name=outSchAssess
#SBATCH --output=./SlurmOutput/outSchAssess-%j.out
#
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --exclude=node[120-128]


julia outageScheduleAssessment.jl ./data/IEEE-RTS96 ./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData_Schedule_MS1-10_Days1-182_v2.csv 16 1 1 1 32 1 182

#julia outageScheduleAssessment.jl ./data/IEEE-RTS96 ./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData_PandzicSchedule_T182.csv 16 1 1 49 64 1 182

#julia outageScheduleAssessment.jl ./data/IEEE-RTS96 ./data/IEEE-RTS96/out_sch_cases/Case1OutageRequestData_LeiWuSchedule_MS1-10_Days1-182_10iters.csv 16 1 1 33 48 1 182
