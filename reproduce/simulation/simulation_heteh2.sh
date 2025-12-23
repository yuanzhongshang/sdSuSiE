#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --job-name=simulation_heteh2
#SBATCH --partition=mulan,main
#SBATCH --exclude=r6406,r6407,r6324
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-100%100
#SBATCH --output=/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/hetehe/out/out_%a.out
#SBATCH --error=/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsexsimulation/hetehe/out/out_%a.err

bash

let k=0
for ((rep=1; rep<=100; rep++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript simulation_heteh2.R ${rep}
fi
done
