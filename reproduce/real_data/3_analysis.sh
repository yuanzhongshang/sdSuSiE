#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --job-name=4079_irnt
#SBATCH --partition=main,mulan
#SBATCH --mem=100G
#SBATCH --array=1-1%200
#SBATCH --output=/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/4079_irnt/out/out_%a.out
#SBATCH --error=/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/4079_irnt/out/out_%a.err

bash

let k=0
for ((rep=1; rep<=17239; rep++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript 3_susie.R ${rep} 4079_irnt
fi
done
